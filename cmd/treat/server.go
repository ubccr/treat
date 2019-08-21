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
	"encoding/gob"
	"fmt"
	"html/template"
	"net/http"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/aebruno/nwalgo"
	"github.com/carbocation/interpose"
	"github.com/gorilla/mux"
	"github.com/gorilla/schema"
	"github.com/gorilla/sessions"
	"github.com/sirupsen/logrus"
	"github.com/ubccr/treat"
)

const (
	TREAT_COOKIE_SESSION = "treat-session"
	TREAT_COOKIE_DB      = "dbname"
	TREAT_COOKIE_SEARCH  = "search"
)

type Application struct {
	templates   map[string]*template.Template
	tmpldir     string
	enableCache bool
	dbs         map[string]*Database
	defaultDb   string
	decoder     *schema.Decoder
	cookieStore *sessions.CookieStore
}

type Database struct {
	name                string
	storage             *Storage
	geneTemplates       map[string]*treat.Template
	geneSamples         map[string][]string
	geneKnockDowns      map[string][]string
	geneReplicates      map[string][]int
	maxEditStop         map[string]int
	maxJuncLen          map[string]int
	maxJuncEnd          map[string]int
	genes               []string
	defaultGene         string
	cache               map[string][]byte
	cacheEditStopTotals map[string]map[int]map[string]float64
}

func init() {
	gob.Register(&SearchFields{})
}

func NewApplication(dbpath, tmpldir string, enableCache bool) (*Application, error) {
	app := &Application{}
	app.dbs = make(map[string]*Database)

	fi, err := os.Stat(dbpath)
	if err != nil {
		return nil, err
	}

	if fi.IsDir() {
		dbabs, _ := filepath.Abs(dbpath)
		dbfiles, err := filepath.Glob(filepath.Join(dbabs, "*.db"))
		if err != nil {
			return nil, err
		}
		for _, d := range dbfiles {
			abs, _ := filepath.Abs(d)
			base := filepath.Base(abs)

			err = app.loadDb(base, abs)
			if err != nil {
				return nil, err
			}

			if len(app.defaultDb) == 0 {
				app.defaultDb = base
			}
		}
	} else {
		abs, _ := filepath.Abs(dbpath)
		base := filepath.Base(abs)
		err = app.loadDb(base, abs)
		if err != nil {
			return nil, err
		}
		app.defaultDb = base
	}

	logrus.Infof("Default DB is: %s", app.defaultDb)

	if len(app.dbs) == 0 {
		return nil, fmt.Errorf("No db files found")
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
		tmpldir = filepath.Join(dir, "/templates")
	}
	logrus.Printf("Using template dir: %s", tmpldir)

	tmpls, err := filepath.Glob(filepath.Join(tmpldir, "*.html"))
	if err != nil {
		return nil, err
	}

	funcMap := template.FuncMap{
		"increment":   incrementFunc,
		"decrement":   decrementFunc,
		"percent":     percent,
		"round":       roundFunc,
		"juncseq":     juncseqFunc,
		"pctSearch":   pctSearchFunc,
		"pctEditStop": pctEditStopFunc,
		"align":       alignFunc,
	}

	app.templates = make(map[string]*template.Template)
	for _, t := range tmpls {
		base := filepath.Base(t)
		if base != "layout.html" && base != "search-form.html" {
			app.templates[base] = template.Must(template.New("layout").Funcs(funcMap).ParseFiles(t,
				filepath.Join(tmpldir, "layout.html"),
				filepath.Join(tmpldir, "search-form.html")))
		}
	}

	app.enableCache = enableCache
	app.tmpldir = tmpldir
	app.decoder = schema.NewDecoder()
	// Create secure cookie with in-secure key. Make this configurable in the future
	app.cookieStore = sessions.NewCookieStore([]byte("not-secure"))
	app.decoder.IgnoreUnknownKeys(true)

	return app, nil
}

func (a *Application) loadDb(base, dbpath string) error {
	logrus.Infof("Processing database: %s", base)
	db := &Database{name: base}
	stg, err := NewStorage(dbpath)
	if err != nil {
		return err
	}
	db.storage = stg

	db.geneTemplates, err = db.storage.TemplateMap()
	if err != nil {
		return err
	}
	if len(db.geneTemplates) == 0 {
		return fmt.Errorf("No genes/templates found. Please load some data first")
	}

	db.cacheEditStopTotals = make(map[string]map[int]map[string]float64)
	db.maxEditStop = make(map[string]int)
	db.maxJuncLen = make(map[string]int)
	db.maxJuncEnd = make(map[string]int)
	db.geneSamples = make(map[string][]string)
	db.geneKnockDowns = make(map[string][]string)
	db.geneReplicates = make(map[string][]int)
	db.genes = make([]string, 0)
	for k := range db.geneTemplates {
		db.genes = append(db.genes, k)

		s, err := db.storage.Samples(k)
		if err != nil {
			return err
		}

		if len(s) == 0 {
			return fmt.Errorf("No samples found for gene %s. Please load some data first", k)
		}

		db.geneSamples[k] = s

		db.geneKnockDowns[k], err = db.storage.KnockDowns(k)
		if err != nil {
			return err
		}

		db.geneReplicates[k], err = db.storage.Replicates(k)
		if err != nil {
			return err
		}

		logrus.Printf("Computing cache for gene %s...", k)
		if _, ok := db.cacheEditStopTotals[k]; !ok {
			db.cacheEditStopTotals[k] = make(map[int]map[string]float64)
		}

		fields := &SearchFields{Gene: k, EditStop: -1, JuncEnd: -1, JuncLen: -1}
		err = db.storage.Search(fields, func(key *treat.AlignmentKey, aln *treat.Alignment) {
			if _, ok := db.cacheEditStopTotals[k][aln.EditStop]; !ok {
				db.cacheEditStopTotals[k][aln.EditStop] = make(map[string]float64)
			}
			db.cacheEditStopTotals[k][aln.EditStop][key.Sample] += aln.Norm

			if aln.EditStop > db.maxEditStop[k] {
				db.maxEditStop[k] = aln.EditStop
			}
			if aln.JuncLen > db.maxJuncLen[k] {
				db.maxJuncLen[k] = aln.JuncLen
			}
			if aln.JuncEnd > db.maxJuncEnd[k] {
				db.maxJuncEnd[k] = aln.JuncEnd
			}
		})

		logrus.Infof("Max ESS: %d Max JL: %d Max JE: %d", db.maxEditStop[k], db.maxJuncLen[k], db.maxJuncEnd[k])
		if db.maxEditStop[k] <= 0 {
			logrus.Errorf("Max edit stop is 0 for gene %s. Expect histogram charts to be broken", k)
		}
		if db.maxJuncLen[k] <= 0 {
			logrus.Errorf("Max junc len is 0 for gene %s. Expect histogram charts to be broken", k)
		}
		if db.maxJuncEnd[k] <= 0 {
			logrus.Errorf("Max junc end is 0 for gene %s. Expect histogram charts to be broken", k)
		}

		if err != nil {
			return fmt.Errorf("Failed computing edit stop totals for gene: %s", k)
		}
	}

	sort.Strings(db.genes)
	// Set default gene for dropdown menu
	db.defaultGene = db.genes[0]

	db.cache = make(map[string][]byte)
	a.dbs[base] = db

	return nil
}

func (a *Application) GetDb(name string) (*Database, error) {
	db, ok := a.dbs[name]
	if !ok {
		return nil, fmt.Errorf("Database not found: %s", name)
	}

	return db, nil
}

func (a *Application) GetDbFromContext(r *http.Request) (*Database, error) {
	dbp := r.Context().Value("db")
	if dbp == nil {
		return nil, fmt.Errorf("Db not found in context")
	}

	db, ok := dbp.(*Database)
	if !ok {
		return nil, fmt.Errorf("Invalid Db pointer in context")
	}

	if db == nil {
		return nil, fmt.Errorf("Db was nil")
	}

	return db, nil
}

func (a *Application) NewSearchFields(w http.ResponseWriter, r *http.Request, db *Database) (*SearchFields, error) {
	vals := r.URL.Query()
	fields := new(SearchFields)
	// set defaults
	fields.EditStop = -2
	fields.JuncLen = -2
	fields.JuncEnd = -2
	fields.Gene = db.defaultGene
	fields.Limit = 10

	session, _ := a.cookieStore.Get(r, TREAT_COOKIE_SESSION)
	search := session.Values[TREAT_COOKIE_SEARCH]
	if search != nil {
		var ok bool
		fields, ok = search.(*SearchFields)
		if !ok {
			fields := new(SearchFields)
			// set defaults
			fields.EditStop = -2
			fields.JuncLen = -2
			fields.JuncEnd = -2
			fields.Gene = db.defaultGene
			fields.Limit = 10
		}
	}

	// Always default to close
	fields.FormOpen = false

	// URL overrides any cookie values
	err := a.decoder.Decode(fields, vals)

	if err != nil {
		return nil, err
	}

	if fields.FormOpen {
		if vals.Get("gene") == "" {
			fields.Gene = db.defaultGene
		}
		if vals.Get("sample") == "" {
			fields.Sample = []string{}
		}
		if vals.Get("kd") == "" {
			fields.KnockDown = []string{}
		}
		if vals.Get("rep") == "" {
			fields.Replicate = []int{}
		}
		if vals.Get("limit") == "" {
			fields.Limit = 10
		}
		if vals.Get("has_mutation") != "1" {
			fields.HasMutation = false
		}
		if vals.Get("has_alt") != "1" {
			fields.HasAlt = false
		}
		if vals.Get("alt") == "" {
			fields.AltRegion = 0
		}
		if vals.Get("edit_stop") == "" {
			fields.EditStop = -2
		}
		if vals.Get("junc_len") == "" {
			fields.JuncLen = -2
		}
		if vals.Get("junc_end") == "" {
			fields.JuncEnd = -2
		}
		if vals.Get("tet") == "" {
			fields.Tetracycline = ""
		}
	}

	session.Values[TREAT_COOKIE_SEARCH] = fields
	err = session.Save(r, w)
	if err != nil {
		logrus.WithFields(logrus.Fields{
			"error": err.Error(),
		}).Error("Failed to update search in session")
	}

	return fields, nil
}

func (a *Application) middlewareStruct() (*interpose.Middleware, error) {
	mw := interpose.New()
	mw.Use(func(next http.Handler) http.Handler { return DbContext(a, next) })
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
	router.Path("/data/bubble").Handler(BubbleJson(a)).Methods("GET")
	router.Path("/data/tmpl").Handler(TemplateSummaryHistogramHandler(a)).Methods("GET")
	router.Path("/heat").Handler(HeatHandler(a)).Methods("GET")
	router.Path("/bubble").Handler(BubbleHandler(a)).Methods("GET")
	router.Path("/search").Handler(SearchHandler(a)).Methods("GET")
	router.Path("/show").Handler(ShowHandler(a)).Methods("GET")
	router.Path("/stats").Handler(StatsHandler(a)).Methods("GET")
	router.Path("/db").Handler(DbHandler(a)).Methods("GET")
	router.Path("/tmpl-report").Handler(TemplateSummaryHandler(a)).Methods("GET")

	return router
}

func writeBase(buf []string, ai int, base rune, count, max uint32, cat string) {
	buf[ai] += `<td class="tcell ` + cat + `">`
	buf[ai] += strings.Repeat("-", int(max-count))
	if count > 0 {
		buf[ai] += strings.Repeat(string(base), int(count))
	}
	buf[ai] += `</td>`
}

func alignFunc(a *treat.Alignment, frag *treat.Fragment, tmpl *treat.Template) template.HTML {

	labels := []string{"FE", "PE"}
	for i := range tmpl.AltRegion {
		labels = append(labels, fmt.Sprintf("A%d", i+1))
	}
	labels = append(labels, "RD")

	aln1, aln2, _ := nwalgo.Align(tmpl.Bases, frag.Bases, 1, -1, -1)

	fragCount := tmpl.Size() + 2
	n := len(aln1)

	buf := make([][]string, fragCount)
	for i := range buf {
		buf[i] = make([]string, n+1)
	}

	fi := 0
	ti := 0
	for ai := 0; ai < n; ai++ {
		hilite := ""
		if n-ti+int(tmpl.EditOffset) == a.EditStop {
			hilite = "hilite"
		}

		if aln1[ai] == '-' {
			buf[0][ai] = `<td class="text-center base-index"></td><td class="text-center base-index"></td>`
			for i := range tmpl.EditSite {
				writeBase(buf[i+1], ai, tmpl.EditBase, 0, frag.EditSite[fi], "ME")
				buf[i+1][ai] += `<td class="text-center base">-</td>`
			}

			writeBase(buf[fragCount-1], ai, frag.EditBase, frag.EditSite[fi], frag.EditSite[fi], "mutant")
			buf[fragCount-1][ai] += `<td class="text-center mutant base">` + string(frag.Bases[fi]) + `</td>`
			fi++
		} else if aln2[ai] == '-' {
			buf[0][ai] = `<td class="text-center ` + hilite + `">` + fmt.Sprintf("%d", tmpl.IndexLabel(n-ti)) + `</td><td class="text-center base-index">` + fmt.Sprintf("%d", tmpl.BaseIndex[ti]) + `</td>`
			max := tmpl.Max(ti)

			for i, t := range tmpl.EditSite {
				writeBase(buf[i+1], ai, tmpl.EditBase, t[ti], max, labels[i])
				buf[i+1][ai] += `<td class="text-center base">` + string(tmpl.Bases[ti]) + `</td>`
			}
			writeBase(buf[fragCount-1], ai, '-', 0, max, "ME")
			buf[fragCount-1][ai] += `<td class="text-center mutant base">-</td>`
			ti++
		} else {
			buf[0][ai] = `<td class="text-center ` + hilite + `">` + fmt.Sprintf("%d", tmpl.IndexLabel(n-ti)) + `</td><td class="text-center base-index">` + fmt.Sprintf("%d", tmpl.BaseIndex[ti]) + `</td>`
			max := tmpl.Max(ti)
			if frag.EditSite[fi] > max {
				max = frag.EditSite[fi]
			}
			cat := ""
			boldi := -1

			if a.AltEditing > 0 && n-ti > tmpl.AltRegion[a.AltEditing-1].Start && n-ti < tmpl.AltRegion[a.AltEditing-1].End {
				boldi = int(a.AltEditing) + 1
				cat = fmt.Sprintf("A%d", a.AltEditing)
			} else if n-ti+int(tmpl.EditOffset) > a.EditStop {
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
				for i, t := range tmpl.EditSite[2:] {
					if frag.EditSite[fi] == t[ti] {
						boldi = i + 2
						cat = fmt.Sprintf("A%d", i+1)
					}
				}
			}

			if n-ti+int(tmpl.EditOffset) > a.EditStop && n-ti+int(tmpl.EditOffset) <= a.JuncEnd {
				cat += " junction"
			}

			for i, t := range tmpl.EditSite {
				bold := ""
				if boldi == i {
					bold = "hilite"
				}
				writeBase(buf[i+1], ai, tmpl.EditBase, t[ti], max, labels[i]+" "+bold)
				buf[i+1][ai] += `<td class="text-center base">` + string(tmpl.Bases[ti]) + `</td>`
			}
			writeBase(buf[fragCount-1], ai, frag.EditBase, frag.EditSite[fi], max, cat)
			buf[fragCount-1][ai] += `<td class="text-center base">` + string(frag.Bases[fi]) + `</td>`
			fi++
			ti++
		}
	}

	// Last edit site has only EditBases
	buf[0][n] = `<td class="text-center">` + fmt.Sprintf("%d", tmpl.IndexLabel(0)) + `</td>`
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
	rows := len(buf[0]) / cols
	if (len(buf[0]) % cols) > 0 {
		rows++
	}

	html := ""

	for r := 0; r < rows; r++ {
		for i, b := range buf {
			html += `<tr>`
			if i == 0 {
				html += `<td>&nbsp;</td>`
			} else {
				html += `<td class="` + labels[i-1] + `">` + labels[i-1] + `</td>`
			}

			end := (r * cols) + cols
			fill := 0
			if end > len(b) {
				fill = end - len(b)
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
		html += `<tr style="border: 0"><td style="border: 0" colspan="` + fmt.Sprintf("%d", (cols*2)+1) + `"></td></tr>`
	}

	return template.HTML(html)
}

func incrementFunc(x int) int {
	x++
	return x
}

func decrementFunc(x int) int {
	x--
	return x
}

func roundFunc(val float64) string {
	return fmt.Sprintf("%.2f", val)
}

func pctSearchFunc(a *treat.Alignment, totals map[string]float64) string {
	y := totals[a.Key.Sample]
	if y == 0 {
		return "0.0"
	}

	d := (a.Norm / y) * 100
	return fmt.Sprintf("%.4f", d)
}

func pctEditStopFunc(a *treat.Alignment, totals map[int]map[string]float64) string {
	y := totals[a.EditStop][a.Key.Sample]
	if y == 0 {
		return "0.0"
	}

	d := (a.Norm / y) * 100
	return fmt.Sprintf("%.4f", d)
}

func juncseqFunc(val string) template.HTML {
	html := ""
	for _, b := range val {
		if b == 'T' || b == 't' {
			html += `<span style="color: red">` + string(b) + `</span>`
		} else {
			html += string(b)
		}
	}

	return template.HTML(html)
}

func Server(dbpath, tmpldir string, port int, enableCache bool) {

	app, err := NewApplication(dbpath, tmpldir, enableCache)
	if err != nil {
		logrus.Fatal(err.Error())
	}

	middle, err := app.middlewareStruct()
	if err != nil {
		logrus.Fatal(err.Error())
	}

	if enableCache {
		logrus.Info("URL caching enabled")
	}

	http.Handle("/", middle)
	logrus.Printf("Running on http://127.0.0.1:%d", port)

	http.ListenAndServe(fmt.Sprintf(":%d", port), nil)
}
