// Copyright 2015 TREAT Authors. All rights reserved.
// Use of this source code is governed by a BSD style
// license that can be found in the LICENSE file.

package main

import (
	"fmt"
	"os"
	"sort"
	"strings"

	"github.com/Sirupsen/logrus"
	"github.com/aebruno/gofasta"
	"github.com/aebruno/nwalgo"
	"github.com/ubccr/treat"
)

// From: https://groups.google.com/d/msg/golang-nuts/FT7cjmcL7gw/Gj4_aEsE_IsJ
type Pair struct {
	Key   string
	Value int
}

type PairList []Pair

func (p PairList) Swap(i, j int)      { p[i], p[j] = p[j], p[i] }
func (p PairList) Len() int           { return len(p) }
func (p PairList) Less(i, j int) bool { return p[i].Value > p[j].Value }

func sortMapByValue(m map[string]int) PairList {
	p := make(PairList, len(m))
	i := 0
	for k, v := range m {
		p[i] = Pair{k, v}
		i++
	}
	sort.Sort(p)
	return p
}

func Mutant(options *AlignOptions, fragments []string, n int) {
	if len(options.TemplatePath) == 0 {
		logrus.Fatal("Please provide path to templates file")
	}
	if len(fragments) == 0 {
		logrus.Fatal("Please provide path to fragment file")
	}
	if len(options.EditBase) != 1 {
		logrus.Fatal("Please provide the edit base")
	}

	tmpl, err := treat.NewTemplateFromFasta(options.TemplatePath, treat.FORWARD, rune(options.EditBase[0]))
	if err != nil {
		logrus.Fatal(err)
	}

	if len(options.Primer3) > 0 {
		err = tmpl.SetPrimer3(options.Primer3)
		if err != nil {
			logrus.Fatal(err)
		}
	}

	if len(options.Primer5) > 0 {
		err = tmpl.SetPrimer5(options.Primer5)
		if err != nil {
			logrus.Fatal(err)
		}
	}

	tm := make(map[string]int)
	fm := make(map[string]int)
	for _, path := range fragments {
		f, err := os.Open(path)
		if err != nil {
			logrus.Fatal(err)
		}
		defer f.Close()

		for rec := range gofasta.SimpleParser(f) {
			frag := treat.NewFragment(rec.Id, rec.Seq, treat.FORWARD, 0, 0, rune(options.EditBase[0]))
			aln1, aln2, _ := nwalgo.Align(tmpl.Bases, frag.Bases, 1, -1, -1)
			if strings.Index(aln1, "-") != -1 {
				tm[aln1]++
			}
			if strings.Index(aln2, "-") != -1 {
				fm[aln2]++
			}
		}
	}

	fmt.Println("Template Indels")
	count := 0
	for _, p := range sortMapByValue(tm) {
		fmt.Printf("%d: %s\n", p.Value, p.Key)
		count++
		if count > n {
			break
		}
	}

	fmt.Println("\nFragment Indels")
	count = 0
	for _, p := range sortMapByValue(fm) {
		fmt.Printf("%d: %s\n", p.Value, p.Key)
		count++
		if count > n {
			break
		}
	}
}
