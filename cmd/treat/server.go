// Copyright 2015 TREAT Authors. All rights reserved.
// Use of this source code is governed by a BSD style
// license that can be found in the LICENSE file.

package main

import (
    "fmt"
    "log"
    "sort"
    "os"
    "bytes"
    "strings"
    "reflect"
    "strconv"
    "encoding/json"
    "encoding/csv"
    "html/template"
    "path/filepath"
    "net/http"
    "net/url"
    "github.com/aebruno/nwalgo"
    "github.com/gorilla/schema"
    "github.com/gorilla/mux"
    "github.com/ubccr/treat"
)


var templates map[string]*template.Template
var db *Storage
var decoder *schema.Decoder
var geneTemplates map[string]*treat.Template
var geneSamples map[string][]string
var maxEditStop map[string]uint32
var maxJuncLen map[string]uint32
var maxJuncEnd map[string]uint32
var genes []string
var cache map[string][]byte
var cacheEditStopTotals map[string]map[uint32]map[string]float64

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

func pctSearch(a *treat.Alignment, totals map[string]float64) (string) {
    y := totals[a.Key.Sample]
    if y == 0 {
        return "0.0"
    }

    d := (a.Norm / y)*100
    return fmt.Sprintf("%.4f", d)
}

func pctEditStop(a *treat.Alignment, totals map[uint32]map[string]float64) (string) {
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

func writeBase(buf []string, ai int, base rune, count, max treat.BaseCountType, cat string) {
    buf[ai] += `<td class="tcell `+cat+`">`
    buf[ai] += strings.Repeat("-", int(max-count))
    if count > 0 {
        buf[ai] += strings.Repeat(string(base), int(count))
    }
    buf[ai] += `</td>`
}

func align(a *treat.Alignment, frag *treat.Fragment, tmpl *treat.Template) (template.HTML) {

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
    if err != nil || frag == nil {
        log.Printf("fragment not found")
        notFoundHandler(w, r)
        return
    }

    alignment, err := db.GetAlignment(key, uint64(id))
    if err != nil || alignment == nil {
        log.Printf("alignment not found")
        notFoundHandler(w, r)
        return
    }

    vars := map[string]interface{}{
        "Template": tmpl,
        "Fragment": frag,
        "Alignment": alignment,
        "Key": key}

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

    totalMap := make(map[string]float64)

    limit := fields.Limit
    fields.Limit = 0
    fields.Offset = 0
    alignments := make([]*treat.Alignment, 0)
    err = db.Search(fields, func (key *treat.AlignmentKey, a *treat.Alignment) {
        a.Key = key
        alignments = append(alignments, a)

        totalMap[a.Key.Sample] += a.Norm
    })

    if err != nil {
        log.Printf("Error fetching alignments for gene: %s", fields.Gene)
        errorHandler(w, r, http.StatusInternalServerError, "")
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
                pctSearch(a, totalMap),
                pctEditStop(a, cacheEditStopTotals[fields.Gene]),
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
        "EditStopTotals": cacheEditStopTotals[fields.Gene],
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

func HeatHandler(w http.ResponseWriter, r *http.Request) {
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

    vars := map[string]interface{}{
        "Template": tmpl,
        "Fields": fields,
        "Samples": geneSamples[fields.Gene],
        "Pages": []int{10,50,100,1000},
        "Genes": genes}

    renderTemplate(w, "heat.html", vars)
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

func HeatMapJson(w http.ResponseWriter, r *http.Request) {
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

    tmpl, ok := geneTemplates[fields.Gene]

    if !ok {
        log.Printf("Error fetching template for gene: %s", fields.Gene)
        errorHandler(w, r, http.StatusInternalServerError, "Gene not found")
        return
    }

    n := tmpl.Len()

    heat := make([][]float64, n)
    for i := 0; i < n; i++ {
        heat[i] = make([]float64, n)
    }

    err = db.Search(fields, func (key *treat.AlignmentKey, a *treat.Alignment) {
        heat[int(a.EditStop)][int(a.JuncLen)] += a.Norm
    })

    if err != nil {
        log.Printf("Fatal error: %s", err)
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
        log.Printf("Error encoding data as json: %s", err)
        http.Error(w, "Fatal system error", http.StatusInternalServerError)
        return
    }

    w.Write(out)
}

func highChartHist(w http.ResponseWriter, r *http.Request, maxMap map[string]uint32, f func(a *treat.Alignment) uint32) {
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

    tmpl, ok := geneTemplates[fields.Gene]
    if !ok {
        log.Printf("Invalid gene: %s", fields.Gene)
        http.Error(w, "Invalid gene", http.StatusInternalServerError)
        return
    }

    samples := make(map[string]map[uint32]float64)

    max := maxMap[fields.Gene]
    err = db.Search(fields, func (key *treat.AlignmentKey, a *treat.Alignment) {
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
        log.Printf("Fatal error: %s", err)
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
        log.Printf("Error encoding data as json: %s", err)
        http.Error(w, "Fatal system error", http.StatusInternalServerError)
        return
    }

    cache[r.URL.String()] = out
    w.Write(out)
    //json.NewEncoder(w).Encode(data)
}

func JuncLenHistogramHandler(w http.ResponseWriter, r *http.Request) {
    highChartHist(w, r, maxJuncLen, func(a *treat.Alignment) uint32 {
        return a.JuncLen
    })
}

func EditHistogramHandler(w http.ResponseWriter, r *http.Request) {
    highChartHist(w, r, maxEditStop, func(a *treat.Alignment) uint32 {
        return a.EditStop
    })
}

func JuncEndHistogramHandler(w http.ResponseWriter, r *http.Request) {
    highChartHist(w, r, maxJuncEnd, func(a *treat.Alignment) uint32 {
        return a.JuncEnd
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

    cacheEditStopTotals = make(map[string]map[uint32]map[string]float64)
    maxEditStop = make(map[string]uint32)
    maxJuncLen = make(map[string]uint32)
    maxJuncEnd = make(map[string]uint32)
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

        log.Printf("Computing edit stop site cache for gene %s...", k)
        if _, ok := cacheEditStopTotals[k]; !ok {
            cacheEditStopTotals[k] = make(map[uint32]map[string]float64)
        }

        fields := &SearchFields{Gene: k, EditStop: -1, JuncEnd: -1, JuncLen: -1}
        err = db.Search(fields, func (key *treat.AlignmentKey, a *treat.Alignment) {
            if _, ok := cacheEditStopTotals[k][a.EditStop]; !ok {
                cacheEditStopTotals[k][a.EditStop] = make(map[string]float64)
            }
            cacheEditStopTotals[k][a.EditStop][key.Sample] += a.Norm

            if a.EditStop > maxEditStop[k] {
                maxEditStop[k] = a.EditStop
            }
            if a.JuncLen > maxJuncLen[k] {
                maxJuncLen[k] = a.JuncLen
            }
            if a.JuncEnd > maxJuncEnd[k] {
                maxJuncEnd[k] = a.JuncEnd
            }
        })

        if err != nil {
            log.Fatalf("Failed computing edit stop totals for gene: %s", k)
        }
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
        "pctSearch": pctSearch,
        "pctEditStop": pctEditStop,
        "grna": grna,
        "align": align,
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

    log.Printf("Running on http://127.0.0.1:%d", port)

    mx := mux.NewRouter()
    mx.NotFoundHandler = http.HandlerFunc(notFoundHandler)
    mx.PathPrefix("/static/").Handler(http.StripPrefix("/static/", http.FileServer(http.Dir(fmt.Sprintf("%s/static", tmpldir)))))
    mx.HandleFunc("/", IndexHandler)
    mx.HandleFunc("/heat", HeatHandler)
    mx.HandleFunc("/search", SearchHandler)
    mx.HandleFunc("/show", ShowHandler)
    mx.HandleFunc("/data/es-hist", EditHistogramHandler)
    mx.HandleFunc("/data/jl-hist", JuncLenHistogramHandler)
    mx.HandleFunc("/data/je-hist", JuncEndHistogramHandler)
    mx.HandleFunc("/data/heat", HeatMapJson)
    http.Handle("/", mx)
    http.ListenAndServe(fmt.Sprintf(":%d", port), nil)
}
