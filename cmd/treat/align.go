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
	"bufio"
	"fmt"
	"os"

	"github.com/Sirupsen/logrus"
	"github.com/aebruno/gofasta"
	"github.com/ubccr/treat"
)

type AlignOptions struct {
	TemplatePath string
	FragmentPath string
	Primer5      string
	Primer3      string
	EditBase     string
	S1           string
	S2           string
}

func PrintAlignment(a1, a2 string, tw int) {
	n := len(a1)
	cols := tw
	rows := n / cols
	if n%cols > 0 {
		rows++
	}

	for r := 0; r < rows; r++ {
		end := (r * cols) + cols
		if end > n {
			end = n
		}

		fmt.Println(a1[(r * cols):end])
		fmt.Println(a2[(r * cols):end])
		fmt.Println()
	}
}

func Align(options *AlignOptions) {
	if len(options.EditBase) != 1 {
		logrus.Fatal("Please provide the edit base")
	}

	if (len(options.S1) > 0 && len(options.S2) == 0) || (len(options.S1) == 0 && len(options.S2) > 0) {
		logrus.Fatal("Please provide 2 fragments to align")
	}
	if len(options.S1) > 0 && len(options.S2) > 0 {
		frag1 := treat.NewFragment("", options.S1, treat.FORWARD, 0, 0, rune(options.EditBase[0]))
		frag2 := treat.NewFragment("", options.S2, treat.FORWARD, 0, 0, rune(options.EditBase[0]))
		aln := new(treat.Alignment)
		a1, a2 := aln.SimpleAlign(frag1, frag2)
		PrintAlignment(a1, a2, 80)
		return
	} else if len(options.FragmentPath) == 0 {
		logrus.Fatal("Please provide either 2 sequences to align or a path to fragment FASTA file")
	}

	var tmpl *treat.Template

	if len(options.TemplatePath) > 0 {
		t, err := treat.NewTemplateFromFasta(options.TemplatePath, treat.FORWARD, rune(options.EditBase[0]))
		if err != nil {
			logrus.Fatal(err)
		}
		tmpl = t

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
	}

	f, err := os.Open(options.FragmentPath)
	if err != nil {
		logrus.Fatal(err)
	}
	defer f.Close()

	if tmpl == nil {
		frags := make([]*treat.Fragment, 0)
		for rec := range gofasta.SimpleParser(f) {
			frag := treat.NewFragment(rec.Id, rec.Seq, treat.FORWARD, 0, 0, rune(options.EditBase[0]))
			frags = append(frags, frag)
			if len(frags) >= 2 {
				break
			}
		}

		if len(frags) < 2 {
			logrus.Fatal("Need 2 fragments to align")
		}

		aln := new(treat.Alignment)
		a1, a2 := aln.SimpleAlign(frags[0], frags[1])
		PrintAlignment(a1, a2, 80)

	} else {
		for rec := range gofasta.SimpleParser(f) {
			frag := treat.NewFragment(rec.Id, rec.Seq, treat.FORWARD, 0, 0, rune(options.EditBase[0]))
			aln := treat.NewAlignment(frag, tmpl, false)
			buf := bufio.NewWriter(os.Stdout)
			aln.WriteTo(buf, frag, tmpl, 80)
			buf.Flush()
		}
	}
}
