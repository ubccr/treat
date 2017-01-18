// Copyright 2015 TREAT Authors. All rights reserved.
//
// This file is part of TREAT.
//
// TREAT is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// TREAT is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with TREAT.  If not, see <http://www.gnu.org/licenses/>.

package main

import (
	"bytes"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"net/http"
	"reflect"
	"sort"
	"strconv"

	"github.com/Sirupsen/logrus"
	"github.com/gorilla/context"
	"github.com/ubccr/treat"
)

func renderTemplate(app *Application, tmpl string, w http.ResponseWriter, data interface{}) {
	if data == nil {
		data = map[string]interface{}{
			"dbs":   app.dbs,
			"curdb": ""}
	}

	var buf bytes.Buffer
	t := app.templates[tmpl]
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
	renderTemplate(app, "error.html", w, nil)
}

func IndexHandler(app *Application) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		db := context.Get(r, "db").(*Database)
		if db == nil {
			logrus.Error("index handler: database not found in request context")
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		fields, err := app.NewSearchFields(w, r, db)

		if err != nil {
			logrus.Printf("Error parsing get request: %s", err)
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		tmpl, ok := db.geneTemplates[fields.Gene]

		if !ok {
			logrus.Printf("Error fetching template for gene: %s", fields.Gene)
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		limit := fields.Limit
		fields.Limit = 0

		count := 0
		err = db.storage.Search(fields, func(key *treat.AlignmentKey, a *treat.Alignment) {
			count++
		})

		fields.Limit = limit

		if err != nil {
			logrus.Printf("Error fetching alignment count for gene: %s", fields.Gene)
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		vars := map[string]interface{}{
			"dbs":      app.dbs,
			"curdb":    db.name,
			"Template": tmpl,
			"Count":    count,
			"Fields":   fields,
			"Samples":  db.geneSamples[fields.Gene],
			"Pages":    []int{10, 50, 100, 1000},
			"Genes":    db.genes}

		renderTemplate(app, "index.html", w, vars)
	})
}

func JuncLenHistogramHandler(app *Application) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		db := context.Get(r, "db").(*Database)
		if db == nil {
			logrus.Error("juncLenHist handler: database not found in request context")
			http.Error(w, "Fatal system error", http.StatusInternalServerError)
			return
		}

		highChartHist(app, w, r, db.maxJuncLen, true, func(a *treat.Alignment) int {
			return a.JuncLen
		})
	})
}

func EditHistogramHandler(app *Application) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		db := context.Get(r, "db").(*Database)
		if db == nil {
			logrus.Error("editHist handler: database not found in request context")
			http.Error(w, "Fatal system error", http.StatusInternalServerError)
			return
		}

		highChartHist(app, w, r, db.maxJuncEnd, false, func(a *treat.Alignment) int {
			return a.EditStop
		})
	})
}

func JuncEndHistogramHandler(app *Application) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		db := context.Get(r, "db").(*Database)
		if db == nil {
			logrus.Error("JunEndHist handler: database not found in request context")
			http.Error(w, "Fatal system error", http.StatusInternalServerError)
			return
		}

		highChartHist(app, w, r, db.maxJuncEnd, false, func(a *treat.Alignment) int {
			return a.JuncEnd
		})
	})
}

func highChartHist(app *Application, w http.ResponseWriter, r *http.Request, maxMap map[string]int, junclen bool, f func(a *treat.Alignment) int) {
	db := context.Get(r, "db").(*Database)
	if db == nil {
		logrus.Error("highChartHist handler: database not found in request context")
		http.Error(w, "Fatal system error", http.StatusInternalServerError)
		return
	}

	if app.enableCache {
		if _, ok := db.cache[r.URL.String()]; ok {
			w.Write(db.cache[r.URL.String()])
			return
		}
	}

	fields, err := app.NewSearchFields(w, r, db)
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

	tmpl, ok := db.geneTemplates[fields.Gene]
	if !ok {
		logrus.Printf("Invalid gene: %s", fields.Gene)
		http.Error(w, "Invalid gene", http.StatusInternalServerError)
		return
	}

	samples := make(map[string]map[int]float64)

	max := maxMap[fields.Gene]
	err = db.storage.Search(fields, func(key *treat.AlignmentKey, a *treat.Alignment) {
		if a.EditStop == int(tmpl.EditStop) && a.JuncLen == 0 {
			return
		}

		if _, ok := samples[key.Sample]; !ok {
			samples[key.Sample] = make(map[int]float64)
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
	var skeys []string
	for k := range samples {
		skeys = append(skeys, k)
	}
	sort.Strings(skeys)
	offset := 0
	if !junclen {
		offset = int(tmpl.EditOffset)
	}
	for _, k := range skeys {
		v := samples[k]
		x := make([]float64, max+1-offset+1)
		for i := range x {
			if _, ok := v[i+offset]; ok {
				x[max-i-offset] = v[i+offset]
			}
		}

		x[max+1-offset] = v[-1+offset]
		m := make(map[string]interface{})
		m["data"] = x
		m["name"] = k
		m["type"] = "spline"
		series = append(series, m)
	}

	cats := make([]int, max+1-offset+1)
	cats[0] = -1 + offset
	for i := range cats {
		cats[i] = int(max) - i
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

		for _, rec := range series {
			for i := len(cats) - 1; i >= 0; i-- {
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

	if app.enableCache {
		db.cache[r.URL.String()] = out
	}

	w.Write(out)
	//json.NewEncoder(w).Encode(data)
}

func ShowHandler(app *Application) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		db := context.Get(r, "db").(*Database)
		if db == nil {
			logrus.Error("show handler: database not found in request context")
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		gene := r.URL.Query().Get("gene")
		tmpl, ok := db.geneTemplates[gene]

		if !ok {
			logrus.Printf("Error fetching template for gene: %s", gene)
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		sample := ""
		for _, s := range db.geneSamples[gene] {
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

		frag, err := db.storage.GetFragment(key, uint64(id))
		if err != nil || frag == nil {
			logrus.Printf("fragment not found: %s", err)
			w.WriteHeader(http.StatusNotFound)
			renderTemplate(app, "404.html", w, nil)
			return
		}

		alignment, err := db.storage.GetAlignment(key, uint64(id))
		if err != nil || alignment == nil {
			logrus.Printf("alignment not found")
			w.WriteHeader(http.StatusNotFound)
			renderTemplate(app, "404.html", w, nil)
			return
		}

		vars := map[string]interface{}{
			"dbs":       app.dbs,
			"curdb":     db.name,
			"Template":  tmpl,
			"Fragment":  frag,
			"Alignment": alignment,
			"Key":       key}

		renderTemplate(app, "show.html", w, vars)
	})
}

func SearchHandler(app *Application) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		db := context.Get(r, "db").(*Database)
		if db == nil {
			logrus.Error("search handler: database not found in request context")
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		fields, err := app.NewSearchFields(w, r, db)

		if err != nil {
			logrus.Printf("Error parsing get request: %s", err)
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		tmpl, ok := db.geneTemplates[fields.Gene]

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
		err = db.storage.Search(fields, func(key *treat.AlignmentKey, a *treat.Alignment) {
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

			for _, a := range alignments {
				csvout.Write([]string{
					strconv.Itoa(int(a.Id)),
					a.Key.Gene,
					a.Key.Sample,
					strconv.Itoa(int(a.ReadCount)),
					fmt.Sprintf("%.4f", a.Norm),
					pctSearchFunc(a, totalMap),
					pctEditStopFunc(a, db.cacheEditStopTotals[fields.Gene]),
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

		if page > (int(len(alignments)/fields.Limit) + 1) {
			page = (int(len(alignments)/fields.Limit) + 1)
		}

		fields.Offset = (page - 1) * fields.Limit
		end := (fields.Offset + fields.Limit)
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
			"dbs":            app.dbs,
			"curdb":          db.name,
			"Template":       tmpl,
			"Count":          count,
			"SearchTotals":   totalMap,
			"EditStopTotals": db.cacheEditStopTotals[fields.Gene],
			"Showing":        showing,
			"Page":           page,
			"Query":          r.URL.RawQuery,
			"Fields":         fields,
			"Alignments":     alignments[fields.Offset:end],
			"Samples":        db.geneSamples[fields.Gene],
			"Pages":          []int{10, 50, 100, 1000},
			"Genes":          db.genes}

		renderTemplate(app, "search.html", w, vars)
	})
}

func HeatHandler(app *Application) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		db := context.Get(r, "db").(*Database)
		if db == nil {
			logrus.Error("heat handler: database not found in request context")
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		fields, err := app.NewSearchFields(w, r, db)

		if err != nil {
			logrus.Printf("Error parsing get request: %s", err)
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		tmpl, ok := db.geneTemplates[fields.Gene]

		if !ok {
			logrus.Printf("Error fetching template for gene: %s", fields.Gene)
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		vars := map[string]interface{}{
			"dbs":      app.dbs,
			"curdb":    db.name,
			"Template": tmpl,
			"Fields":   fields,
			"Samples":  db.geneSamples[fields.Gene],
			"Pages":    []int{10, 50, 100, 1000},
			"Genes":    db.genes}

		renderTemplate(app, "heat.html", w, vars)
	})
}

func HeatMapJson(app *Application) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		db := context.Get(r, "db").(*Database)
		if db == nil {
			logrus.Error("heatmap json handler: database not found in request context")
			http.Error(w, "Fatal error", http.StatusInternalServerError)
			return
		}

		fields, err := app.NewSearchFields(w, r, db)
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

		tmpl, ok := db.geneTemplates[fields.Gene]

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

		err = db.storage.Search(fields, func(key *treat.AlignmentKey, a *treat.Alignment) {
			if a.EditStop >= int(tmpl.EditOffset) {
				heat[a.EditStop-int(tmpl.EditOffset)][a.JuncLen] += a.Norm
			}
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
				series[k][0] = i + int(tmpl.EditOffset)
				series[k][1] = j
				if i+int(tmpl.EditOffset) == int(tmpl.EditStop) && j == 0 {
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

func BubbleHandler(app *Application) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		db := context.Get(r, "db").(*Database)
		if db == nil {
			logrus.Error("bubble handler: database not found in request context")
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		fields, err := app.NewSearchFields(w, r, db)

		if err != nil {
			logrus.Printf("Error parsing get request: %s", err)
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		tmpl, ok := db.geneTemplates[fields.Gene]

		if !ok {
			logrus.Printf("Error fetching template for gene: %s", fields.Gene)
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		vars := map[string]interface{}{
			"dbs":      app.dbs,
			"curdb":    db.name,
			"Template": tmpl,
			"Fields":   fields,
			"Samples":  db.geneSamples[fields.Gene],
			"Pages":    []int{10, 50, 100, 1000},
			"Genes":    db.genes}

		renderTemplate(app, "bubble.html", w, vars)
	})
}

func BubbleJson(app *Application) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		db := context.Get(r, "db").(*Database)
		if db == nil {
			logrus.Error("bubble json handler: database not found in request context")
			http.Error(w, "Fatal error", http.StatusInternalServerError)
			return
		}

		fields, err := app.NewSearchFields(w, r, db)
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

		tmpl, ok := db.geneTemplates[fields.Gene]

		if !ok {
			logrus.Printf("Error fetching template for gene: %s", fields.Gene)
			http.Error(w, "Gene not found", http.StatusInternalServerError)
			return
		}

		bubbleMap := make(map[int]map[uint32]float64)

		err = db.storage.Search(fields, func(key *treat.AlignmentKey, a *treat.Alignment) {
			if a.EditStop >= int(tmpl.EditOffset) {
				frag, err := db.storage.GetFragment(key, a.Id)
				if err != nil || frag == nil {
					logrus.Printf("fragment not found: %s", err)
					return
				}

				for i, t := range frag.EditSite {
					eb, ok := bubbleMap[i]
					if !ok {
						eb = make(map[uint32]float64)
					}
					eb[t] += a.Norm
					bubbleMap[i] = eb
				}
			}
		})

		if err != nil {
			logrus.Printf("Fatal error: %s", err)
			http.Error(w, "Fatal database error.", http.StatusInternalServerError)
			return
		}

		type editSiteBubble struct {
			Name  int             `json:"name"`
			Total float64         `json:"total"`
			T     [][]interface{} `json:"t"`
			Pre   int             `json:"pre"`
			Full  int             `json:"full"`
		}

		n := tmpl.Len()
		bubbles := make([]*editSiteBubble, n)

		for key, val := range bubbleMap {
			b := &editSiteBubble{
				Name: (n - 1) - int(key) + int(tmpl.EditOffset),
				T:    make([][]interface{}, 0),
				Full: int(tmpl.EditSite[0][key]),
				Pre:  int(tmpl.EditSite[1][key]),
			}

			for t, cnt := range val {
				b.T = append(b.T, []interface{}{int(t), cnt})
				b.Total += cnt
			}
			bubbles[int(key)] = b
		}

		out, err := json.Marshal(bubbles)
		if err != nil {
			logrus.Printf("Error encoding bubble data as json: %s", err)
			http.Error(w, "Fatal system error", http.StatusInternalServerError)
			return
		}

		w.Write(out)
	})
}

func TemplateSummaryHandler(app *Application) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		db := context.Get(r, "db").(*Database)
		if db == nil {
			logrus.Error("tmpl summary handler: database not found in request context")
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		fields, err := app.NewSearchFields(w, r, db)

		if err != nil {
			logrus.Printf("Error parsing get request: %s", err)
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		tmpl, ok := db.geneTemplates[fields.Gene]

		if !ok {
			logrus.Printf("Error fetching template for gene: %s", fields.Gene)
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		vars := map[string]interface{}{
			"dbs":      app.dbs,
			"curdb":    db.name,
			"Template": tmpl,
			"Fields":   fields,
			"Samples":  db.geneSamples[fields.Gene],
			"Pages":    []int{10, 50, 100, 1000},
			"Genes":    db.genes}

		renderTemplate(app, "tmpl-report.html", w, vars)
	})
}

func TemplateSummaryHistogramHandler(app *Application) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		db := context.Get(r, "db").(*Database)
		if db == nil {
			logrus.Error("tmpl summary json handler: database not found in request context")
			http.Error(w, "Fatal error", http.StatusInternalServerError)
			return
		}

		fields, err := app.NewSearchFields(w, r, db)
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

		tmpl, ok := db.geneTemplates[fields.Gene]
		if !ok {
			logrus.Printf("Invalid gene: %s", fields.Gene)
			http.Error(w, "Invalid gene", http.StatusInternalServerError)
			return
		}

		samples := make(map[string]map[int]float64)

		err = db.storage.Search(fields, func(key *treat.AlignmentKey, a *treat.Alignment) {
			ok := false
			if a.EditStop == int(tmpl.EditStop) && a.JuncLen == 0 {
				ok = true
			} else if a.EditStop == tmpl.Len()-1+int(tmpl.EditOffset) && a.JuncLen == 0 {
				ok = true
			}

			if !ok {
				return
			}

			if _, ok = samples[key.Sample]; !ok {
				samples[key.Sample] = make(map[int]float64)
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

		var skeys []string
		for k := range samples {
			skeys = append(skeys, k)
		}
		sort.Strings(skeys)
		for _, k := range skeys {
			v := samples[k]
			x := []float64{v[int(tmpl.EditStop)]}
			y := []float64{v[tmpl.Len()-1+int(tmpl.EditOffset)]}

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

func DbHandler(app *Application) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		session, _ := app.cookieStore.Get(r, TREAT_COOKIE_SESSION)

		name := r.URL.Query().Get("name")
		_, err := app.GetDb(name)
		if err != nil {
			logrus.WithFields(logrus.Fields{
				"dbname": name,
			}).Error("Invalid database name")
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		session.Values[TREAT_COOKIE_DB] = name
		err = session.Save(r, w)
		if err != nil {
			logrus.WithFields(logrus.Fields{
				"error": err.Error(),
			}).Error("Db handler: failed to save session")
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		http.Redirect(w, r, "/", 302)
	})
}

func StatsHandler(app *Application) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		db := context.Get(r, "db").(*Database)
		if db == nil {
			logrus.Error("stats handler: database not found in request context")
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		countByString := r.FormValue("countby")
		countBy := COUNT_FRAG
		if "Unique" == countByString {
			countBy = COUNT_UNIQUE
		} else {
			countByString = "Fragments"
		}

		fields, err := app.NewSearchFields(w, r, db)
		if err != nil {
			logrus.Printf("Error parsing get request: %s", err)
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		tmpl, ok := db.geneTemplates[fields.Gene]
		if !ok {
			logrus.Printf("Error fetching template for gene: %s", fields.Gene)
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		stats, err := geneStats(db.storage, fields.Gene, countBy)
		if err != nil {
			logrus.Printf("Failed to compute stats for gene %s: %s", fields.Gene, err)
			errorHandler(app, w, http.StatusInternalServerError)
			return
		}

		vars := map[string]interface{}{
			"dbs":      app.dbs,
			"curdb":    db.name,
			"stats":    stats,
			"Fields":   fields,
			"Template": tmpl,
			"Counts":   []string{"Fragments", "Unique"},
			"Countby":  countByString,
			"Genes":    db.genes}

		renderTemplate(app, "stats.html", w, vars)
	})
}
