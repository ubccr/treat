package main

import (
    "fmt"
    "log"
    "os"
    "encoding/json"
    "html/template"
    "path/filepath"
    "net/http"
    "github.com/gorilla/mux"
    "github.com/ubccr/treat"
)

var templates map[string]*template.Template
var db *treat.Storage

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
    aln := make([]*treat.Alignment, 0)
    fields := &treat.SearchFields{Limit: 10, All: true}

    err := db.Search(fields, func (k string, a *treat.Alignment) {
        aln = append(aln, a)
    })

    if err != nil {
        log.Printf("Fatal error: %s", err)
        http.Error(w, "Fatal database error.", http.StatusInternalServerError)
        return
    }

    json.NewEncoder(w).Encode(aln)
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
