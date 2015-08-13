// Copyright 2015 TREAT Authors. All rights reserved.
// Use of this source code is governed by a BSD style
// license that can be found in the LICENSE file.

package main

import (
    "bytes"
    "html/template"
    "net/http"
    "encoding/json"
    "encoding/csv"
    "reflect"
    "strconv"
    "fmt"
    "sort"

    "github.com/Sirupsen/logrus"
    "github.com/ubccr/treat"
)


func renderTemplate(w http.ResponseWriter, t *template.Template, data interface{}) {
    var buf bytes.Buffer
    err := t.ExecuteTemplate(&buf, "layout", data)

    if err != nil {
        logrus.Printf("Error rendering template: %s", err)
        http.Error(w, "Fatal error rendering template", http.StatusInternalServerError)
        return
    }

    buf.WriteTo(w)
}

func errorHandler(app *Application, w http.ResponseWriter, status int) {
    w.WriteHeader(status)

    renderTemplate(w, app.templates["error.html"], nil)
}

func IndexHandler(app *Application) http.Handler {
    return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
        fields, err := app.NewSearchFields(r.URL)

        if err != nil {
            logrus.Printf("Error parsing get request: %s", err)
            errorHandler(app, w, http.StatusInternalServerError)
            return
        }

        tmpl,ok := app.geneTemplates[fields.Gene]

        if !ok {
            logrus.Printf("Error fetching template for gene: %s", fields.Gene)
            errorHandler(app, w, http.StatusInternalServerError)
            return
        }

        count := 0
        err = app.db.Search(fields, func (key *treat.AlignmentKey, a *treat.Alignment) {
            count++
        })

        if err != nil {
            logrus.Printf("Error fetching alignment count for gene: %s", fields.Gene)
            errorHandler(app, w, http.StatusInternalServerError)
            return
        }

        vars := map[string]interface{}{
            "Template": tmpl,
            "Count": count,
            "Fields": fields,
            "Samples": app.geneSamples[fields.Gene],
            "Pages": []int{10,50,100,1000},
            "Genes": app.genes}

        renderTemplate(w, app.templates["index.html"], vars)
    })
}

func JuncLenHistogramHandler(app *Application) http.Handler {
    return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
        highChartHist(app, w, r, app.maxJuncEnd, func(a *treat.Alignment) uint32 {
            return a.JuncLen
        })
    })
}

func EditHistogramHandler(app *Application) http.Handler {
    return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
        highChartHist(app, w, r, app.maxJuncEnd, func(a *treat.Alignment) uint32 {
            return a.EditStop
        })
    })
}

func JuncEndHistogramHandler(app *Application) http.Handler {
    return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
        highChartHist(app, w, r, app.maxJuncEnd, func(a *treat.Alignment) uint32 {
            return a.JuncEnd
        })
    })
}

func highChartHist(app *Application, w http.ResponseWriter, r *http.Request, maxMap map[string]uint32, f func(a *treat.Alignment) uint32) {
    if _, ok := app.cache[r.URL.String()]; ok {
        w.Write(app.cache[r.URL.String()])
        return
    }

    fields, err := app.NewSearchFields(r.URL)
    fields.Limit = 0
    fields.Offset = 0

    if err != nil {
        logrus.Printf("Error parsing get request: %s", err)
        http.Error(w, "Invalid get parameter in request", http.StatusInternalServerError)
        return
    }

    if len(fields.Gene) == 0 {
        logrus.Printf("Missing required field gene")
        http.Error(w, "Missing required field gene", http.StatusInternalServerError)
        return
    }

    tmpl, ok := app.geneTemplates[fields.Gene]
    if !ok {
        logrus.Printf("Invalid gene: %s", fields.Gene)
        http.Error(w, "Invalid gene", http.StatusInternalServerError)
        return
    }

    samples := make(map[string]map[uint32]float64)

    max := maxMap[fields.Gene]
    err = app.db.Search(fields, func (key *treat.AlignmentKey, a *treat.Alignment) {
        if a.EditStop == tmpl.EditStop && a.JuncLen == 0 {
            return
        }

        if _, ok := samples[key.Sample]; !ok {
            samples[key.Sample] = make(map[uint32]float64)
        }

        val := f(a)

        samples[key.Sample][val] += a.Norm
    })

    if err != nil {
        logrus.Printf("Fatal error: %s", err)
        http.Error(w, "Fatal database error.", http.StatusInternalServerError)
        return
    }

    series := make([]map[string]interface{}, 0)
    for k,v := range(samples) {
        x := make([]float64, max+1)
        for i := range(x) {
            if _, ok := v[uint32(i)]; ok {
                x[int(max)-i] = v[uint32(i)]
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
        defer csvout.Flush()
        col := "edit_stop"
        if r.URL.Path == "/data/jl-hist" {
            col = "junc_len"
        } else if r.URL.Path == "/data/je-hist" {
            col = "junc_end"
        }
        csvout.Write([]string{col, "name", "norm_count"})

        w.Header().Set("Content-Type", "text/csv; charset=utf-8")
        w.Header().Set("Content-Disposition", "attachment; filename="+col+".csv")

        for _, rec := range(series) {
            for i := len(cats)-1; i >= 0; i-- {
                es := cats[i]
                norm := reflect.ValueOf(rec["data"])
                name := reflect.ValueOf(rec["name"])
                csvout.Write([]string{strconv.Itoa(es), fmt.Sprintf("%s", name), fmt.Sprintf("%.4f", norm.Index(i).Float())})
            }
        }

        return
    }

    out, err := json.Marshal(data)
    if err != nil {
        logrus.Printf("Error encoding data as json: %s", err)
        http.Error(w, "Fatal system error", http.StatusInternalServerError)
        return
    }

    app.cache[r.URL.String()] = out
    w.Write(out)
    //json.NewEncoder(w).Encode(data)
}

func ShowHandler(app *Application) http.Handler {
    return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
        gene := r.URL.Query().Get("gene")
        tmpl,ok := app.geneTemplates[gene]

        if !ok {
            logrus.Printf("Error fetching template for gene: %s", gene)
            errorHandler(app, w, http.StatusInternalServerError)
            return
        }

        sample := ""
        for _, s := range(app.geneSamples[gene]) {
            if s == r.URL.Query().Get("sample") {
                sample = s
            }
        }

        if len(sample) == 0 {
            logrus.Printf("sample not found: %s", r.URL.Query().Get("sample"))
            errorHandler(app, w, http.StatusInternalServerError)
            return
        }

        id, err := strconv.Atoi(r.URL.Query().Get("id"))
        if err != nil {
            logrus.Printf("id not found: %s", r.URL.Query().Get("id"))
            errorHandler(app, w, http.StatusInternalServerError)
            return
        }

        key := &treat.AlignmentKey{Gene: gene, Sample: sample}

        frag, err := app.db.GetFragment(key, uint64(id))
        if err != nil || frag == nil {
            logrus.Printf("fragment not found")
            w.WriteHeader(http.StatusNotFound)
            renderTemplate(w, app.templates["404.html"], nil)
            return
        }

        alignment, err := app.db.GetAlignment(key, uint64(id))
        if err != nil || alignment == nil {
            logrus.Printf("alignment not found")
            w.WriteHeader(http.StatusNotFound)
            renderTemplate(w, app.templates["404.html"], nil)
            return
        }

        vars := map[string]interface{}{
            "Template": tmpl,
            "Fragment": frag,
            "Alignment": alignment,
            "Key": key}

        renderTemplate(w, app.templates["show.html"], vars)
    })
}

func SearchHandler(app *Application) http.Handler {
    return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
        fields, err := app.NewSearchFields(r.URL)

        if err != nil {
            logrus.Printf("Error parsing get request: %s", err)
            errorHandler(app, w, http.StatusInternalServerError)
            return
        }

        tmpl,ok := app.geneTemplates[fields.Gene]

        if !ok {
            logrus.Printf("Error fetching template for gene: %s", fields.Gene)
            errorHandler(app, w, http.StatusInternalServerError)
            return
        }

        totalMap := make(map[string]float64)

        limit := fields.Limit
        fields.Limit = 0
        fields.Offset = 0
        alignments := make([]*treat.Alignment, 0)
        err = app.db.Search(fields, func (key *treat.AlignmentKey, a *treat.Alignment) {
            a.Key = key
            alignments = append(alignments, a)

            totalMap[a.Key.Sample] += a.Norm
        })

        if err != nil {
            logrus.Printf("Error fetching alignments for gene: %s", fields.Gene)
            errorHandler(app, w, http.StatusInternalServerError)
            return
        }

        sort.Sort(ByReadCount{alignments})

        if r.URL.Query().Get("export") == "1" {
            csvout := csv.NewWriter(w)
            defer csvout.Flush()
            csvout.Write([]string{"id", "gene", "sample", "read_count", "norm_count", "pct_search", "pct_edit_stop", "edit_stop", "junc_end", "junc_len", "junc_seq"})

            w.Header().Set("Content-Type", "text/csv; charset=utf-8")
            w.Header().Set("Content-Disposition", "attachment; filename=treat-export.csv")

            for _, a := range(alignments) {
                csvout.Write([]string{
                    strconv.Itoa(int(a.Id)),
                    a.Key.Gene,
                    a.Key.Sample,
                    strconv.Itoa(int(a.ReadCount)),
                    fmt.Sprintf("%.4f", a.Norm),
                    pctSearchFunc(a, totalMap),
                    pctEditStopFunc(a, app.cacheEditStopTotals[fields.Gene]),
                    strconv.Itoa(int(a.EditStop)),
                    strconv.Itoa(int(a.JuncEnd)),
                    strconv.Itoa(int(a.JuncLen)),
                    a.JuncSeq})
            }

            return
        }

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
        if end > len(alignments) {
            end = len(alignments)
        }
        if end < 0 {
            end = 0
        }

        count := len(alignments)
        showing := fields.Offset + fields.Limit
        if showing > count {
            showing = count
        }

        vars := map[string]interface{}{
            "Template": tmpl,
            "Count": count,
            "SearchTotals": totalMap,
            "EditStopTotals": app.cacheEditStopTotals[fields.Gene],
            "Showing": showing,
            "Page": page,
            "Query": r.URL.RawQuery,
            "Fields": fields,
            "Alignments": alignments[fields.Offset:end],
            "Samples": app.geneSamples[fields.Gene],
            "Pages": []int{10,50,100,1000},
            "Genes": app.genes}

        renderTemplate(w, app.templates["search.html"], vars)
    })
}

func HeatHandler(app *Application) http.Handler {
    return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
        fields, err := app.NewSearchFields(r.URL)

        if err != nil {
            logrus.Printf("Error parsing get request: %s", err)
            errorHandler(app, w, http.StatusInternalServerError)
            return
        }

        tmpl,ok := app.geneTemplates[fields.Gene]

        if !ok {
            logrus.Printf("Error fetching template for gene: %s", fields.Gene)
            errorHandler(app, w, http.StatusInternalServerError)
            return
        }

        vars := map[string]interface{}{
            "Template": tmpl,
            "Fields": fields,
            "Samples": app.geneSamples[fields.Gene],
            "Pages": []int{10,50,100,1000},
            "Genes": app.genes}

        renderTemplate(w, app.templates["heat.html"], vars)
    })
}

func HeatMapJson(app *Application) http.Handler {
    return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
        fields, err := app.NewSearchFields(r.URL)
        fields.Limit = 0
        fields.Offset = 0

        if err != nil {
            logrus.Printf("Error parsing get request: %s", err)
            http.Error(w, "Invalid get parameter in request", http.StatusInternalServerError)
            return
        }

        if len(fields.Gene) == 0 {
            logrus.Printf("Missing required field gene")
            http.Error(w, "Missing required field gene", http.StatusInternalServerError)
            return
        }

        tmpl, ok := app.geneTemplates[fields.Gene]

        if !ok {
            logrus.Printf("Error fetching template for gene: %s", fields.Gene)
            http.Error(w, "Gene not found", http.StatusInternalServerError)
            return
        }

        n := tmpl.Len()

        heat := make([][]float64, n)
        for i := 0; i < n; i++ {
            heat[i] = make([]float64, n)
        }

        err = app.db.Search(fields, func (key *treat.AlignmentKey, a *treat.Alignment) {
            heat[int(a.EditStop)][int(a.JuncLen)] += a.Norm
        })

        if err != nil {
            logrus.Printf("Fatal error: %s", err)
            http.Error(w, "Fatal database error.", http.StatusInternalServerError)
            return
        }

        series := make([][]interface{}, n*n)
        max := float64(0.0)
        k := 0
        for i := 0; i < n; i++ {
            for j := 0; j < n; j++ {
                series[k] = make([]interface{}, 3)
                series[k][0] = i
                series[k][1] = j
                if i == int(tmpl.EditStop) && j == 0 {
                    series[k][2] = 0.0
                } else {
                    series[k][2] = heat[i][j]
                    if heat[i][j] > max {
                        max = heat[i][j]
                    }
                }
                k++
            }
        }

        data := make(map[string]interface{})
        data["series"] = series
        data["max"] = max

        out, err := json.Marshal(data)
        if err != nil {
            logrus.Printf("Error encoding data as json: %s", err)
            http.Error(w, "Fatal system error", http.StatusInternalServerError)
            return
        }

        w.Write(out)
    })
}

func TemplateSummaryHandler(app *Application) http.Handler {
    return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
        fields, err := app.NewSearchFields(r.URL)

        if err != nil {
            logrus.Printf("Error parsing get request: %s", err)
            errorHandler(app, w, http.StatusInternalServerError)
            return
        }

        tmpl,ok := app.geneTemplates[fields.Gene]

        if !ok {
            logrus.Printf("Error fetching template for gene: %s", fields.Gene)
            errorHandler(app, w, http.StatusInternalServerError)
            return
        }

        vars := map[string]interface{}{
            "Template": tmpl,
            "Fields": fields,
            "Samples": app.geneSamples[fields.Gene],
            "Pages": []int{10,50,100,1000},
            "Genes": app.genes}

        renderTemplate(w, app.templates["tmpl-report.html"], vars)
    })
}

func TemplateSummaryHistogramHandler(app *Application) http.Handler {
    return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
        fields, err := app.NewSearchFields(r.URL)
        fields.Limit = 0
        fields.Offset = 0

        if err != nil {
            logrus.Printf("Error parsing get request: %s", err)
            http.Error(w, "Invalid get parameter in request", http.StatusInternalServerError)
            return
        }

        if len(fields.Gene) == 0 {
            logrus.Printf("Missing required field gene")
            http.Error(w, "Missing required field gene", http.StatusInternalServerError)
            return
        }

        tmpl, ok := app.geneTemplates[fields.Gene]
        if !ok {
            logrus.Printf("Invalid gene: %s", fields.Gene)
            http.Error(w, "Invalid gene", http.StatusInternalServerError)
            return
        }

        samples := make(map[string]map[uint32]float64)
        primer5 := (tmpl.Len()-2)-tmpl.Primer5

        err = app.db.Search(fields, func (key *treat.AlignmentKey, a *treat.Alignment) {
            ok := false
            if a.EditStop == tmpl.EditStop && a.JuncLen == 0 {
                ok = true
            } else if a.EditStop == uint32(primer5) && a.JuncLen == 0 {
                ok = true
            }

            if !ok {
                return
            }

            if _, ok = samples[key.Sample]; !ok {
                samples[key.Sample] = make(map[uint32]float64)
            }

            samples[key.Sample][a.EditStop] += a.Norm
        })

        if err != nil {
            logrus.Printf("Fatal error: %s", err)
            http.Error(w, "Fatal database error.", http.StatusInternalServerError)
            return
        }

        fe := make([]map[string]interface{}, 0)
        pe := make([]map[string]interface{}, 0)
        for k,v := range(samples) {
            x := []float64{v[tmpl.EditStop]}
            y := []float64{v[uint32(primer5)]}

            m := make(map[string]interface{})
            m["data"] = y
            m["name"] = k
            fe = append(fe, m)

            m = make(map[string]interface{})
            m["data"] = x
            m["name"] = k
            pe = append(pe, m)
        }

        data := make(map[string]interface{})
        data["pe"] = pe
        data["fe"] = fe

        out, err := json.Marshal(data)
        if err != nil {
            logrus.Printf("Error encoding data as json: %s", err)
            http.Error(w, "Fatal system error", http.StatusInternalServerError)
            return
        }

        w.Write(out)
    })
}
