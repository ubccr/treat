package main

import (
    "fmt"
    "log"
    "os"
    "strings"
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

func renderTemplate(w http.ResponseWriter, name string, data interface{}) {
    err := templates[name].ExecuteTemplate(w, "layout", data)

    if err != nil {
        http.Error(w, err.Error(), http.StatusInternalServerError)
    }
}

func IndexHandler(w http.ResponseWriter, r *http.Request) {
    renderTemplate(w, "index.html", nil)
}

func EditHistogramHandler(w http.ResponseWriter, r *http.Request) {
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

    err = db.Search(fields, func (k string, a *treat.Alignment) {
        key := strings.Split(string(k), ";")

        if _, ok := samples[key[1]]; !ok {
            samples[key[1]] = make([]float64, tmpl.Len())
        }

        samples[key[1]][a.EditStop] += a.Norm
    })

    if err != nil {
        log.Printf("Fatal error: %s", err)
        http.Error(w, "Fatal database error.", http.StatusInternalServerError)
        return
    }

    json.NewEncoder(w).Encode(samples)
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
        if base != "layout.html" {
            templates[base] = template.Must(template.New("layout").ParseFiles(t, tmpldir + "/templates/layout.html"))
        }
    }

    mx := mux.NewRouter()
    mx.PathPrefix("/static/").Handler(http.StripPrefix("/static/", http.FileServer(http.Dir(fmt.Sprintf("%s/static", tmpldir)))))
    mx.HandleFunc("/", IndexHandler)
    mx.HandleFunc("/data/es-hist", EditHistogramHandler)
    http.Handle("/", mx)
    http.ListenAndServe(fmt.Sprintf(":%d", port), nil)
}
