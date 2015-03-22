package main

import (
    "fmt"
    "log"
    "sort"
    "os"
    "bytes"
    "reflect"
    "strconv"
    "encoding/json"
    "encoding/csv"
    "html/template"
    "path/filepath"
    "net/http"
    "net/url"
    "github.com/gorilla/schema"
    "github.com/gorilla/mux"
    "github.com/ubccr/treat"
)


var templates map[string]*template.Template
var db *Storage
var decoder *schema.Decoder
var geneTemplates map[string]*treat.Template
var geneSamples map[string][]string
var genes []string
var cache map[string][]byte

func increment(x int) (int) {
    x++
    return x
}

func decrement(x int) (int) {
    x--
    return x
}

func round(val float64) (string) {
    return fmt.Sprintf("%.4f", val)
}

func pct(a *treat.Alignment, totals map[uint64]map[string]float64) (string) {
    y := totals[a.EditStop][a.Key.Sample]
    if y == 0 {
        return "0.0"
    }

    d := (a.Norm / y)*100
    return fmt.Sprintf("%.4f", d)
}

func grna(a *treat.Alignment) (template.HTML) {
    html := ""
    for i := 0; i < a.GrnaEdit.BitLen(); i++ {
        if a.GrnaEdit.Bit(i) == 1 {
            html += `<span class="label label-success">`+fmt.Sprintf("gRNA%d", i+1)+`</span> `
        }
    }
    return template.HTML(html)
}

func juncseq(val string) (template.HTML) {
    html := ""
    for _,b := range(val) {
        if b == 'T' || b == 't' {
            html += `<span style="color: red">`+string(b)+`</span>`
        } else {
            html += string(b)
        }
    }

    return template.HTML(html)
}

func renderTemplate(w http.ResponseWriter, name string, data interface{}) {
    var buf bytes.Buffer
    err := templates[name].ExecuteTemplate(&buf, "layout", data)

    if err != nil {
        log.Printf("Error rendering template: %s", err)
        http.Error(w, "Fatal error rendering template", http.StatusInternalServerError)
        return
    }

    buf.WriteTo(w)
}

func notFoundHandler(w http.ResponseWriter, r *http.Request) {
    w.WriteHeader(http.StatusNotFound)
    renderTemplate(w, "404.html", nil)
}

func errorHandler(w http.ResponseWriter, r *http.Request, status int, message string) {
    w.WriteHeader(status)
    renderTemplate(w, "error.html", message)
}

func ShowHandler(w http.ResponseWriter, r *http.Request) {
    gene := r.URL.Query().Get("gene")
    tmpl,ok := geneTemplates[gene]

    if !ok {
        log.Printf("Error fetching template for gene: %s", gene)
        errorHandler(w, r, http.StatusInternalServerError, "Gene not found")
        return
    }

    sample := ""
    for _, s := range(geneSamples[gene]) {
        if s == r.URL.Query().Get("sample") {
            sample = s
        }
    }

    if len(sample) == 0 {
        log.Printf("sample not found: %s", r.URL.Query().Get("sample"))
        errorHandler(w, r, http.StatusInternalServerError, "Sample not found")
        return
    }

    id, err := strconv.Atoi(r.URL.Query().Get("id"))
    if err != nil {
        log.Printf("id not found: %s", r.URL.Query().Get("id"))
        errorHandler(w, r, http.StatusInternalServerError, "Id not found")
        return
    }

    key := &treat.AlignmentKey{Gene: gene, Sample: sample}

    frag, err := db.GetFragment(key, uint64(id))
    if err != nil {
        log.Printf("fragment not found")
        notFoundHandler(w, r)
        return
    }

    vars := map[string]interface{}{
        "Template": tmpl,
        "Gene": gene,
        "Fragment": frag,
        "Sample": sample}

    renderTemplate(w, "show.html", vars)
}

func SearchHandler(w http.ResponseWriter, r *http.Request) {
    fields, err := NewSearchFields(r.URL)

    if err != nil {
        log.Printf("Error parsing get request: %s", err)
        errorHandler(w, r, http.StatusInternalServerError, "")
        return
    }

    tmpl,ok := geneTemplates[fields.Gene]

    if !ok {
        log.Printf("Error fetching template for gene: %s", fields.Gene)
        errorHandler(w, r, http.StatusInternalServerError, "No templates found for gene")
        return
    }

    totalMap := make(map[uint64]map[string]float64)

    limit := fields.Limit
    fields.Limit = 0
    fields.Offset = 0
    alignments := make([]*treat.Alignment, 0)
    err = db.Search(fields, func (key *treat.AlignmentKey, a *treat.Alignment) {
        a.Key = key
        alignments = append(alignments, a)

        if _, ok := totalMap[a.EditStop]; !ok {
            totalMap[a.EditStop] = make(map[string]float64)
        }
        totalMap[a.EditStop][a.Key.Sample] += a.Norm
    })

    if err != nil {
        log.Printf("Error fetching alignments for gene: %s", fields.Gene)
        errorHandler(w, r, http.StatusInternalServerError, "")
        return
    }

    sort.Sort(ByReadCount{alignments})

    fields.Limit = limit
    if fields.Limit == 0 {
        fields.Limit = 10
    }

    page, err := strconv.Atoi(r.URL.Query().Get("page"))
    if err != nil {
        page = 0
    }

    if page <= 0 {
        page = 1
    }

    if page > (int(len(alignments)/fields.Limit)+1) {
        page = (int(len(alignments)/fields.Limit)+1)
    }

    fields.Offset = (page-1)*fields.Limit
    end := (fields.Offset+fields.Limit)
    if end > len(alignments)-1 {
        end = len(alignments)-1
    }
    if end < 0 {
        end = 0
    }

    count := len(alignments)
    showing := fields.Offset + fields.Limit

    vars := map[string]interface{}{
        "Template": tmpl,
        "Count": count,
        "Totals": totalMap,
        "Showing": showing,
        "Page": page,
        "Query": r.URL.RawQuery,
        "Fields": fields,
        "Alignments": alignments[fields.Offset:end],
        "Samples": geneSamples[fields.Gene],
        "Pages": []int{10,50,100,1000},
        "Genes": genes}

    renderTemplate(w, "search.html", vars)
}

func IndexHandler(w http.ResponseWriter, r *http.Request) {
    fields, err := NewSearchFields(r.URL)

    if err != nil {
        log.Printf("Error parsing get request: %s", err)
        errorHandler(w, r, http.StatusInternalServerError, "")
        return
    }

    tmpl,ok := geneTemplates[fields.Gene]

    if !ok {
        log.Printf("Error fetching template for gene: %s", fields.Gene)
        errorHandler(w, r, http.StatusInternalServerError, "No templates found for gene")
        return
    }

    count := 0
    err = db.Search(fields, func (key *treat.AlignmentKey, a *treat.Alignment) {
        count++
    })

    if err != nil {
        log.Printf("Error fetching alignment count for gene: %s", fields.Gene)
        errorHandler(w, r, http.StatusInternalServerError, "")
        return
    }

    vars := map[string]interface{}{
        "Template": tmpl,
        "Count": count,
        "Fields": fields,
        "Samples": geneSamples[fields.Gene],
        "Pages": []int{10,50,100,1000},
        "Genes": genes}

    renderTemplate(w, "index.html", vars)
}

func NewSearchFields(url *url.URL) (*SearchFields, error) {
    vals := url.Query()
    fields := new(SearchFields)
    err := decoder.Decode(fields, vals)

    if err != nil {
        return nil, err
    }

    if len(fields.Gene) == 0 {
        for k := range(geneTemplates) {
            fields.Gene = k
            break
        }
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
    fields.Limit = 0
    fields.Offset = 0

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

    if r.URL.Query().Get("export") == "1" {
        csvout := csv.NewWriter(w)
        col := "edit_stop"
        if r.URL.Path == "/data/jl-hist" {
            col = "junc_len"
        }
        csvout.Write([]string{col, "name", "norm_count"})

        w.Header().Set("Content-Type", "text/csv; charset=utf-8")
        w.Header().Set("Content-Disposition", "attachment; filename="+col+".csv")

        for _, rec := range(series) {
            for i, es := range(cats) {
                norm := reflect.ValueOf(rec["data"])
                name := reflect.ValueOf(rec["name"])
                csvout.Write([]string{strconv.Itoa(es), fmt.Sprintf("%s", name), fmt.Sprintf("%.4f", norm.Index(i).Float())})
            }
        }

        csvout.Flush()
        return
    }

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
    dbx, err := NewStorage(dbpath)
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
        "decrement": decrement,
        "round": round,
        "juncseq": juncseq,
        "pct": pct,
        "grna": grna,
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
    decoder = schema.NewDecoder()
    decoder.IgnoreUnknownKeys(true)

    mx := mux.NewRouter()
    mx.NotFoundHandler = http.HandlerFunc(notFoundHandler)
    mx.PathPrefix("/static/").Handler(http.StripPrefix("/static/", http.FileServer(http.Dir(fmt.Sprintf("%s/static", tmpldir)))))
    mx.HandleFunc("/", IndexHandler)
    mx.HandleFunc("/search", SearchHandler)
    mx.HandleFunc("/show", ShowHandler)
    mx.HandleFunc("/data/es-hist", EditHistogramHandler)
    mx.HandleFunc("/data/jl-hist", JuncHistogramHandler)
    http.Handle("/", mx)
    http.ListenAndServe(fmt.Sprintf(":%d", port), nil)
}
