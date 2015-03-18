package main

import (
    "fmt"
    "log"
    "os"
    "encoding/json"
    "html/template"
    "path/filepath"
    "net/http"
    "github.com/gorilla/schema"
    "github.com/gorilla/mux"
    "github.com/ubccr/treat"
)

var templates map[string]*template.Template
var db *treat.Storage
var decoder = schema.NewDecoder()
var cache map[string][]byte

func renderTemplate(w http.ResponseWriter, name string, data interface{}) {
    err := templates[name].ExecuteTemplate(w, "layout", data)

    if err != nil {
        http.Error(w, err.Error(), http.StatusInternalServerError)
    }
}

func IndexHandler(w http.ResponseWriter, r *http.Request) {
    fields := new(treat.SearchFields)
    err := decoder.Decode(fields, r.URL.Query())

    if err != nil {
        log.Printf("Error parsing get request: %s", err)
        http.Error(w, "Invalid get parameter in request", http.StatusInternalServerError)
        return
    }

    genes, err := db.Genes()
    if err != nil || len(genes) == 0 {
        log.Printf("Error fetching genes: %s", err)
        http.Error(w, "No genes", http.StatusInternalServerError)
        return
    }

    if len(fields.Gene) == 0 {
        fields.Gene = genes[0]
    }

    tmpl, err := db.GetTemplate(fields.Gene)
    if err != nil {
        log.Printf("Error fetching template for gene: %s", fields.Gene)
        http.Error(w, "No templates found for gene", http.StatusInternalServerError)
        return
    }

    vars := map[string]interface{}{
        "Template": tmpl,
        "Genes": genes}

    renderTemplate(w, "index.html", vars)
}

func JuncHistogramHandler(w http.ResponseWriter, r *http.Request) {
    if _, ok := cache[r.URL.String()]; ok {
        w.Write(cache[r.URL.String()])
        return
    }

    fields := new(treat.SearchFields)
    err := decoder.Decode(fields, r.URL.Query())

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

    samples := make(map[string]map[uint64]float64)

    max := uint64(0)
    err = db.Search(fields, func (key *treat.AlignmentKey, a *treat.Alignment) {
        if _, ok := samples[key.Sample]; !ok {
            samples[key.Sample] = make(map[uint64]float64)
        }

        if a.JuncLen > max {
            max = a.JuncLen
        }
        samples[key.Sample][a.JuncLen] += a.Norm
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

func EditHistogramHandler(w http.ResponseWriter, r *http.Request) {
    if _, ok := cache[r.URL.String()]; ok {
        w.Write(cache[r.URL.String()])
        return
    }

    fields := new(treat.SearchFields)
    err := decoder.Decode(fields, r.URL.Query())

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

    tmpl, err := db.GetTemplate(fields.Gene)
    if err != nil {
        log.Printf("Error fetching template for gene: %s", fields.Gene)
        http.Error(w, "No templates found for gene", http.StatusInternalServerError)
        return
    }

    samples := make(map[string][]float64)

    err = db.Search(fields, func (key *treat.AlignmentKey, a *treat.Alignment) {
        if _, ok := samples[key.Sample]; !ok {
            samples[key.Sample] = make([]float64, tmpl.Len())
        }

        samples[key.Sample][a.EditStop] += a.Norm
    })

    if err != nil {
        log.Printf("Fatal error: %s", err)
        http.Error(w, "Fatal database error.", http.StatusInternalServerError)
        return
    }

    series := make([]map[string]interface{}, 0)
    for k,v := range(samples) {
        for i, j := 0, len(v)-1; i < j; i, j = i+1, j-1 {
            v[i], v[j] = v[j], v[i]
        }
        m := make(map[string]interface{})
        m["data"] = v
        m["name"] = k
        m["type"] = "spline"
        series = append(series, m)
    }


    cats := make([]int, tmpl.Len())
    for i := range(cats) {
        cats[i] = (tmpl.Len()-1)-i
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

func Server(dbpath, tmpldir string, port int) {
    dbx, err := treat.NewStorage(dbpath)
    if err != nil {
        log.Fatal(err)
    }
    db = dbx

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

    templates = make(map[string]*template.Template)
    for _, t := range tmpls {
        base := filepath.Base(t)
        if base != "layout.html" && base != "search-form.html" {
            templates[base] = template.Must(template.New("layout").ParseFiles(t,
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
