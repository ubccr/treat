// Copyright 2015 TREAT Authors. All rights reserved.
// Use of this source code is governed by a BSD style
// license that can be found in the LICENSE file.

package main

import (
    "fmt"
    "log"
    "strings"
    "os"
    "html/template"
    "path/filepath"
    "net/http"
    "net/url"

    "github.com/aebruno/nwalgo"
    "github.com/Sirupsen/logrus"
    "github.com/ubccr/treat"
    "github.com/carbocation/interpose"
    "github.com/gorilla/mux"
    "github.com/gorilla/schema"
)

type Application struct {
    templates            map[string]*template.Template
    tmpldir              string
    dbpaths              map[string]string
    curdb                string
    db                   *Storage
    decoder              *schema.Decoder
    geneTemplates        map[string]*treat.Template
    geneSamples          map[string][]string
    maxEditStop          map[string]uint32
    maxJuncLen           map[string]uint32
    maxJuncEnd           map[string]uint32
    genes                []string
    cache                map[string][]byte
    cacheEditStopTotals  map[string]map[uint32]map[string]float64
}

func NewApplication(dbpath, tmpldir string) (*Application, error) {
    app := &Application{}
    app.dbpaths = make(map[string]string)

    fi, err := os.Stat(dbpath)
    if err != nil {
        return nil, err
    }

    if fi.IsDir() {
        dbabs,_ := filepath.Abs(dbpath)
        dbfiles, err := filepath.Glob(dbabs + "/*.db")
        if err != nil {
            return nil, err
        }
        for _, d := range dbfiles {
            abs,_ := filepath.Abs(d)
            base := filepath.Base(abs)
            app.dbpaths[base] = abs

            if len(app.curdb) == 0 {
                app.curdb = base
            }
        }
    } else {
        abs,_ := filepath.Abs(dbpath)
        base := filepath.Base(abs)
        app.dbpaths[base] = abs
        app.curdb = base
    }

    if len(app.dbpaths) == 0 {
        return nil, fmt.Errorf("No db files found")
    }

    err = app.loadDb(app.dbpaths[app.curdb])
    if err != nil {
        return nil, err
    }

    if len(tmpldir) == 0 {
        // default to directory of current executable 
        path, err := filepath.EvalSymlinks(os.Args[0])
        if err != nil {
            return nil, err
        }
        dir, err := filepath.Abs(filepath.Dir(path))
        if err != nil {
            return nil, err
        }
        tmpldir = dir
    }
    log.Printf("Using template dir: %s\n", tmpldir)

    tmpls, err := filepath.Glob(tmpldir + "/templates/*.html")
    if err != nil {
        return nil, err
    }

    funcMap := template.FuncMap{
        "increment": incrementFunc,
        "decrement": decrementFunc,
        "round": roundFunc,
        "juncseq": juncseqFunc,
        "pctSearch": pctSearchFunc,
        "pctEditStop": pctEditStopFunc,
        "align": alignFunc,
    }

    app.templates = make(map[string]*template.Template)
    for _, t := range tmpls {
        base := filepath.Base(t)
        if base != "layout.html" && base != "search-form.html" {
            app.templates[base] = template.Must(template.New("layout").Funcs(funcMap).ParseFiles(t,
                                                        tmpldir + "/templates/layout.html",
                                                        tmpldir + "/templates/search-form.html"))
        }
    }

    app.tmpldir = tmpldir
    app.decoder = schema.NewDecoder()
    app.decoder.IgnoreUnknownKeys(true)

    return app, nil
}

func (a *Application) loadDb(dbpath string) (error) {
    db, err := NewStorage(dbpath)
    if err != nil {
        return err
    }

    a.db = db

    a.geneTemplates, err = a.db.TemplateMap()
    if err != nil {
        return err
    }
    if len(a.geneTemplates) == 0 {
        return fmt.Errorf("No genes/templates found. Please load some data first")
    }

    a.cacheEditStopTotals = make(map[string]map[uint32]map[string]float64)
    a.maxEditStop = make(map[string]uint32)
    a.maxJuncLen = make(map[string]uint32)
    a.maxJuncEnd = make(map[string]uint32)
    a.geneSamples = make(map[string][]string)
    a.genes = make([]string, 0)
    for k := range(a.geneTemplates) {
        a.genes = append(a.genes, k)

        s, err := a.db.Samples(k)
        if err != nil {
            return err
        }

        if len(s) == 0 {
            return fmt.Errorf("No samples found for gene %s. Please load some data first", k)
        }

        a.geneSamples[k] = s

        log.Printf("Computing edit stop site cache for gene %s...", k)
        if _, ok := a.cacheEditStopTotals[k]; !ok {
            a.cacheEditStopTotals[k] = make(map[uint32]map[string]float64)
        }

        fields := &SearchFields{Gene: k, EditStop: -1, JuncEnd: -1, JuncLen: -1}
        err = a.db.Search(fields, func (key *treat.AlignmentKey, aln *treat.Alignment) {
            if _, ok := a.cacheEditStopTotals[k][aln.EditStop]; !ok {
                a.cacheEditStopTotals[k][aln.EditStop] = make(map[string]float64)
            }
            a.cacheEditStopTotals[k][aln.EditStop][key.Sample] += aln.Norm

            if aln.EditStop > a.maxEditStop[k] {
                a.maxEditStop[k] = aln.EditStop
            }
            if aln.JuncLen > a.maxJuncLen[k] {
                a.maxJuncLen[k] = aln.JuncLen
            }
            if aln.JuncEnd > a.maxJuncEnd[k] {
                a.maxJuncEnd[k] = aln.JuncEnd
            }
        })

        if err != nil {
            return fmt.Errorf("Failed computing edit stop totals for gene: %s", k)
        }
    }

    a.cache = make(map[string][]byte)

    return nil
}

func (a *Application) NewSearchFields(url *url.URL) (*SearchFields, error) {
    vals := url.Query()
    fields := new(SearchFields)
    err := a.decoder.Decode(fields, vals)

    if err != nil {
        return nil, err
    }

    if len(fields.Gene) == 0 {
        for k := range(a.geneTemplates) {
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

func (a *Application) middlewareStruct() (*interpose.Middleware, error) {
    mw := interpose.New()
    mw.UseHandler(a.router())

    return mw, nil
}

func (a *Application) router() *mux.Router {
    router := mux.NewRouter()

    router.NotFoundHandler = http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
        w.WriteHeader(http.StatusNotFound)
        renderTemplate(a, "404.html", w, nil)
    })
    router.PathPrefix("/static/").Handler(http.StripPrefix("/static/", http.FileServer(http.Dir(fmt.Sprintf("%s/static", a.tmpldir)))))

    router.Path("/").Handler(IndexHandler(a)).Methods("GET")
    router.Path("/data/es-hist").Handler(EditHistogramHandler(a)).Methods("GET")
    router.Path("/data/jl-hist").Handler(JuncLenHistogramHandler(a)).Methods("GET")
    router.Path("/data/je-hist").Handler(JuncEndHistogramHandler(a)).Methods("GET")
    router.Path("/data/heat").Handler(HeatMapJson(a)).Methods("GET")
    router.Path("/data/tmpl").Handler(TemplateSummaryHistogramHandler(a)).Methods("GET")
    router.Path("/heat").Handler(HeatHandler(a)).Methods("GET")
    router.Path("/search").Handler(SearchHandler(a)).Methods("GET")
    router.Path("/show").Handler(ShowHandler(a)).Methods("GET")
    router.Path("/db").Handler(DbHandler(a)).Methods("GET")
    router.Path("/tmpl-report").Handler(TemplateSummaryHandler(a)).Methods("GET")

    return router
}

func writeBase(buf []string, ai int, base rune, count, max treat.BaseCountType, cat string) {
    buf[ai] += `<td class="tcell `+cat+`">`
    buf[ai] += strings.Repeat("-", int(max-count))
    if count > 0 {
        buf[ai] += strings.Repeat(string(base), int(count))
    }
    buf[ai] += `</td>`
}

func alignFunc(a *treat.Alignment, frag *treat.Fragment, tmpl *treat.Template) (template.HTML) {

    labels := []string{"FE","PE"}
    for i := range(tmpl.AltRegion) {
        labels = append(labels, fmt.Sprintf("A%d", i+1))
    }
    labels = append(labels, "CL")


    aln1, aln2, _ := nwalgo.Align(tmpl.Bases, frag.Bases, 1, -1, -1)

    fragCount := tmpl.Size()+2
    n := len(aln1)

    buf := make([][]string, fragCount)
    for i := range(buf) {
        buf[i] = make([]string, n+1)
    }

    fi := 0
    ti := 0
    for ai := 0; ai < n; ai++ {
        hilite := ""
        if uint32(n-ti) == a.EditStop {
            hilite = "hilite"
        }

        if aln1[ai] == '-' {
            buf[0][ai] = `<td class="text-center base-index"></td><td class="text-center base-index"></td>`
            for i, _ := range tmpl.EditSite {
                writeBase(buf[i+1], ai, tmpl.EditBase, 0, frag.EditSite[fi], "ME")
                buf[i+1][ai] += `<td class="text-center base">-</td>`
            }

            writeBase(buf[fragCount-1], ai, frag.EditBase, frag.EditSite[fi], frag.EditSite[fi], "mutant")
            buf[fragCount-1][ai] += `<td class="text-center mutant base">`+string(frag.Bases[fi])+`</td>`
            fi++
        } else if aln2[ai] == '-' {
            buf[0][ai] = `<td class="text-center `+hilite+`">`+fmt.Sprintf("%d", n-ti)+`</td><td class="text-center base-index">`+fmt.Sprintf("%d", tmpl.BaseIndex[ti])+`</td>`
            max := tmpl.Max(ti)

            for i, t := range tmpl.EditSite {
                writeBase(buf[i+1], ai, tmpl.EditBase, t[ti], max, labels[i])
                buf[i+1][ai] += `<td class="text-center base">`+string(tmpl.Bases[ti])+`</td>`
            }
            writeBase(buf[fragCount-1], ai, '-', 0, max, "ME")
            buf[fragCount-1][ai] += `<td class="text-center mutant base">-</td>`
            ti++
        } else {
            buf[0][ai] = `<td class="text-center `+hilite+`">`+fmt.Sprintf("%d", n-ti)+`</td><td class="text-center base-index">`+fmt.Sprintf("%d", tmpl.BaseIndex[ti])+`</td>`
            max := tmpl.Max(ti)
            if frag.EditSite[fi] > max {
                max = frag.EditSite[fi]
            }
            cat := ""
            boldi := -1
            if uint32(n-ti) > a.EditStop {
                if frag.EditSite[fi] == tmpl.EditSite[1][ti] {
                    cat = "PE"
                    boldi = 1
                } else if frag.EditSite[fi] == tmpl.EditSite[0][ti] {
                    cat = "FE"
                    boldi = 0
                }
            } else {
                if frag.EditSite[fi] == tmpl.EditSite[0][ti] {
                    cat = "FE"
                    boldi = 0
                } else if frag.EditSite[fi] == tmpl.EditSite[1][ti] {
                    cat = "PE"
                    boldi = 1
                }
            }

            if boldi == -1 {
                for i,t := range(tmpl.EditSite[2:]) {
                    if frag.EditSite[fi] == t[ti] {
                        boldi = i+2
                        cat = fmt.Sprintf("A%d", i+1)
                    }
                }
            }

            if uint32(n-ti) > a.EditStop && uint32(n-ti) <= a.JuncEnd {
                cat += " junction"
            }

            for i, t := range tmpl.EditSite {
                bold := ""
                if boldi == i {
                    bold = "hilite"
                }
                writeBase(buf[i+1], ai, tmpl.EditBase, t[ti], max, labels[i]+" "+bold)
                buf[i+1][ai] += `<td class="text-center base">`+string(tmpl.Bases[ti])+`</td>`
            }
            writeBase(buf[fragCount-1], ai, frag.EditBase, frag.EditSite[fi], max, cat)
            buf[fragCount-1][ai] += `<td class="text-center base">`+string(frag.Bases[fi])+`</td>`
            fi++
            ti++
        }
    }

    // Last edit site has only EditBases
    buf[0][n] = `<td class="text-center">`+fmt.Sprintf("%d", 0)+`</td>`
    max := tmpl.Max(ti)
    if frag.EditSite[fi] > max {
        max = frag.EditSite[fi]
    }
    cat := "PE"
    if frag.EditSite[fi] == tmpl.EditSite[0][ti] {
        cat = "FE"
    }

    for i, t := range tmpl.EditSite {
        writeBase(buf[i+1], n, tmpl.EditBase, t[ti], max, labels[i])
    }
    writeBase(buf[fragCount-1], n, frag.EditBase, frag.EditSite[fi], max, cat)

    cols := 17
    rows := len(buf[0])/cols
    if (len(buf[0]) % cols) > 0 {
        rows++
    }

    html := ""

    for r := 0; r < rows; r++ {
        for i, b := range(buf) {
            html += `<tr>`
            if i == 0 {
                html += `<td>&nbsp;</td>`
            } else {
                html += `<td class="`+labels[i-1]+`">`+labels[i-1]+`</td>`
            }

            end := (r*cols)+cols
            fill := 0
            if end > len(b) {
                fill = end-len(b)
                end = len(b)
            }
            html += strings.Join(b[(r*cols):end], "")
            if fill > 0 {
                for x := 0; x < (fill*2)+1; x++ {
                    html += `<td style="border:none">&nbsp</td>`
                }
            }
            html += `</tr>`
        }
        html += `<tr style="border: 0"><td style="border: 0" colspan="`+fmt.Sprintf("%d", (cols*2)+1)+`"></td></tr>`
    }

    return template.HTML(html)
}

func incrementFunc(x int) (int) {
    x++
    return x
}

func decrementFunc(x int) (int) {
    x--
    return x
}

func roundFunc(val float64) (string) {
    return fmt.Sprintf("%.4f", val)
}

func pctSearchFunc(a *treat.Alignment, totals map[string]float64) (string) {
    y := totals[a.Key.Sample]
    if y == 0 {
        return "0.0"
    }

    d := (a.Norm / y)*100
    return fmt.Sprintf("%.4f", d)
}

func pctEditStopFunc(a *treat.Alignment, totals map[uint32]map[string]float64) (string) {
    y := totals[a.EditStop][a.Key.Sample]
    if y == 0 {
        return "0.0"
    }

    d := (a.Norm / y)*100
    return fmt.Sprintf("%.4f", d)
}

func juncseqFunc(val string) (template.HTML) {
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


func Server(dbpath, tmpldir string, port int) {

    app, err := NewApplication(dbpath, tmpldir)
    if err != nil {
        logrus.Fatal(err.Error())
    }

    middle, err := app.middlewareStruct()
    if err != nil {
        logrus.Fatal(err.Error())
    }

    http.Handle("/", middle)
    logrus.Printf("Running on http://127.0.0.1:%d", port)

    http.ListenAndServe(fmt.Sprintf(":%d", port), nil)
}
