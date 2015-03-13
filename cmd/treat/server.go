package main

import (
    "fmt"
    "html/template"
    "net/http"
    "github.com/gorilla/mux"
)

var templates *template.Template

func IndexHandler(w http.ResponseWriter, r *http.Request) {
    templates.ExecuteTemplate(w, "layout", "Index")
}

func Server(tmpldir string, port int) {
    templates = template.Must(template.ParseGlob("./templates/*"))
    mx := mux.NewRouter()
    mx.PathPrefix("/static/").Handler(http.StripPrefix("/static/", http.FileServer(http.Dir("./static"))))
    mx.HandleFunc("/", IndexHandler)
    http.Handle("/", mx)
    http.ListenAndServe(fmt.Sprintf(":%d", port), nil)
}
