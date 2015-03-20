package main

import (
    "fmt"
    "log"
    "os"
    "bytes"
    "encoding/json"
    "html/template"
    "path/filepath"
    "net/http"
    "net/url"
    "github.com/gorilla/schema"
    "github.com/gorilla/mux"
    "github.com/ubccr/treat"
)

var templates map[string]*template.Template
var db *treat.Storage
var decoder = schema.NewDecoder()
var geneTemplates map[string]*treat.Template
var geneSamples map[string][]string
var genes []string
var cache map[string][]byte

func increment(x int) (int) {
    x++
    return x
}

func renderTemplate(w http.ResponseWriter, name string, data interface{}) {
    var buf bytes.Buffer
    err := templates[name].ExecuteTemplate(&buf, "layout", data)

    if err != nil {
        http.Error(w, err.Error(), http.StatusInternalServerError)
        return
    }

    buf.WriteTo(w)
}

func IndexHandler(w http.ResponseWriter, r *http.Request) {
    fields, err := NewSearchFields(r.URL)

    if err != nil {
        log.Printf("Error parsing get request: %s", err)
        http.Error(w, "Invalid get parameter in request", http.StatusInternalServerError)
        return
    }

    if len(fields.Gene) == 0 {
        for k := range(geneTemplates) {
            fields.Gene = k
            break
        }
    }

    tmpl,ok := geneTemplates[fields.Gene]

    if !ok {
        log.Printf("Error fetching template for gene: %s", fields.Gene)
        http.Error(w, "No templates found for gene", http.StatusInternalServerError)
        return
    }

    count := 0
    err = db.Search(fields, func (key *treat.AlignmentKey, a *treat.Alignment) {
        count++
    })

    if err != nil {
        log.Printf("Error fetching alignment count for gene: %s", fields.Gene)
        http.Error(w, "System error", http.StatusInternalServerError)
        return
    }

    vars := map[string]interface{}{
        "Template": tmpl,
        "Count": count,
        "Query": r.URL.RawQuery,
        "Fields": fields,
        "Samples": geneSamples[fields.Gene],
        "Pages": []int{10,50,100,1000},
        "Genes": genes}

    renderTemplate(w, "index.html", vars)
}

func NewSearchFields(url *url.URL) (* treat.SearchFields, error) {
    vals := url.Query()
    fields := new(treat.SearchFields)
    err := decoder.Decode(fields, vals)

    if err != nil {
        return nil, err
    }

    if vals.Get("edit_stop") == "" {
        fields.EditStop = -1
    }
    if vals.Get("junc_len") == "" {
        fields.JuncLen = -1
    }
    if vals.Get("junc_end") == "" {
        fields.JuncEnd = -1
    }

    return fields, nil
}

func highChartHist(w http.ResponseWriter, r *http.Request, f func(a *treat.Alignment) uint64) {
    if _, ok := cache[r.URL.String()]; ok {
        w.Write(cache[r.URL.String()])
        return
    }

    fields, err := NewSearchFields(r.URL)

    if err != nil {
        log.Printf("Error parsing get request: %s", err)
        http.Error(w, "Invalid get parameter in request", http.StatusInternalServerError)
        return
    }

    if len(fields.Gene) == 0 {
        log.Printf("Missing required field gene")
        http.Error(w, "Missing required field gene", http.StatusInternalServerError)
        return
    }

    if _,ok := geneTemplates[fields.Gene]; !ok {
        log.Printf("Invalid gene: %s", fields.Gene)
        http.Error(w, "Invalid gene", http.StatusInternalServerError)
        return
    }

    samples := make(map[string]map[uint64]float64)

    max := uint64(0)
    err = db.Search(fields, func (key *treat.AlignmentKey, a *treat.Alignment) {
        if _, ok := samples[key.Sample]; !ok {
            samples[key.Sample] = make(map[uint64]float64)
        }

        val := f(a)

        if val > max {
            max = val
        }
        samples[key.Sample][val] += a.Norm
    })

    if err != nil {
        log.Printf("Fatal error: %s", err)
        http.Error(w, "Fatal database error.", http.StatusInternalServerError)
        return
    }

    series := make([]map[string]interface{}, 0)
    for k,v := range(samples) {
        x := make([]float64, max+1)
        for i := range(x) {
            if _, ok := v[uint64(i)]; ok {
                x[int(max)-i] = v[uint64(i)]
            }
        }

        m := make(map[string]interface{})
        m["data"] = x
        m["name"] = k
        m["type"] = "spline"
        series = append(series, m)
    }


    cats := make([]int, max+1)
    for i := range(cats) {
        cats[i] = int(max)-i
    }

    data := make(map[string]interface{})
    data["cats"] = cats
    data["series"] = series

    out, err := json.Marshal(data)
    if err != nil {
        log.Printf("Error encoding data as json: %s", err)
        http.Error(w, "Fatal system error", http.StatusInternalServerError)
        return
    }

    cache[r.URL.String()] = out
    w.Write(out)
    //json.NewEncoder(w).Encode(data)
}

func JuncHistogramHandler(w http.ResponseWriter, r *http.Request) {
    highChartHist(w, r, func(a *treat.Alignment) uint64 {
        return a.JuncLen
    })
}

func EditHistogramHandler(w http.ResponseWriter, r *http.Request) {
    highChartHist(w, r, func(a *treat.Alignment) uint64 {
        return a.EditStop
    })
}

func Server(dbpath, tmpldir string, port int) {
    dbx, err := treat.NewStorage(dbpath)
    if err != nil {
        log.Fatal(err)
    }
    db = dbx

    geneTemplates, err = db.TemplateMap()
    if err != nil {
        log.Fatal(err)
    }
    if len(geneTemplates) == 0 {
        log.Fatal("No genes/templates found. Please load some data first")
    }

    geneSamples = make(map[string][]string)
    genes = make([]string, 0)
    for k := range(geneTemplates) {
        genes = append(genes, k)

        s, err := db.Samples(k)
        if err != nil {
            log.Fatal(err)
        }

        if len(s) == 0 {
            log.Fatalf("No samples found for gene %s. Please load some data first", k)
        }

        geneSamples[k] = s
    }

    if len(tmpldir) == 0 {
        // default to directory of current executable 
        path, err := filepath.EvalSymlinks(os.Args[0])
        if err != nil {
            log.Fatal(err)
        }
        dir, err := filepath.Abs(filepath.Dir(path))
        if err != nil {
            log.Fatal(err)
        }
        tmpldir = dir
    }
    log.Printf("Using template dir: %s\n", tmpldir)

    tmpls, err := filepath.Glob(tmpldir + "/templates/*.html")
    if err != nil {
        log.Fatal(err)
    }

    funcMap := template.FuncMap{
        "increment": increment,
    }


    templates = make(map[string]*template.Template)
    for _, t := range tmpls {
        base := filepath.Base(t)
        if base != "layout.html" && base != "search-form.html" {
            templates[base] = template.Must(template.New("layout").Funcs(funcMap).ParseFiles(t,
                                                        tmpldir + "/templates/layout.html",
                                                        tmpldir + "/templates/search-form.html"))
        }
    }

    cache = make(map[string][]byte)

    mx := mux.NewRouter()
    mx.PathPrefix("/static/").Handler(http.StripPrefix("/static/", http.FileServer(http.Dir(fmt.Sprintf("%s/static", tmpldir)))))
    mx.HandleFunc("/", IndexHandler)
    mx.HandleFunc("/data/es-hist", EditHistogramHandler)
    mx.HandleFunc("/data/jl-hist", JuncHistogramHandler)
    http.Handle("/", mx)
    http.ListenAndServe(fmt.Sprintf(":%d", port), nil)
}
