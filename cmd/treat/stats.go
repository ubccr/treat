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
	"fmt"
	"strings"

	"github.com/Sirupsen/logrus"
	"github.com/ubccr/treat"
)

func percent(x, y int) float64 {
	return (float64(x) / float64(y)) * float64(100)
}

func Stats(dbpath, gene string) {
	s, err := NewStorage(dbpath)
	if err != nil {
		logrus.Fatal(err)
	}
	fmt.Printf("db path: %s\n", dbpath)
	fmt.Printf("version: %.1f\n\n", s.version)

	geneTemplates, err := s.TemplateMap()
	if err != nil {
		logrus.Fatal(err)
	}

	for g, tmpl := range geneTemplates {
		if len(gene) > 0 && g != gene {
			continue
		}

		grandTotal := 0
		nonMutant := 0
		countMap := make(map[string]int)
		totalMap := make(map[string]int)
		err = s.Search(&SearchFields{Gene: g, All: true, EditStop: -1, JuncLen: -1, JuncEnd: -1}, func(key *treat.AlignmentKey, a *treat.Alignment) {
			if a.HasMutation == uint8(0) {
				countMap[key.Sample]++
				nonMutant++
			}

			totalMap[key.Sample]++
			grandTotal++
		})

		if err != nil {
			logrus.Fatal(err)
		}

		fmt.Println(strings.Repeat("=", 80))
		fmt.Println(g)
		fmt.Println(strings.Repeat("=", 80))
		fmt.Printf("%20s%11d\n", "Total Alignments:", grandTotal)
		fmt.Printf("%20s%11d\n", "Non-Mutant:", nonMutant)
		fmt.Printf("%20s%11d\n", "Mutant:", grandTotal-nonMutant)
		fmt.Printf("%20s%11d\n", "Edit Stop:", tmpl.EditStop)
		fmt.Printf("%20s%11s\n", "Edit Base:", string(tmpl.EditBase))
		fmt.Printf("%20s%11d\n", "Alt Templates:", len(tmpl.AltRegion))
		fmt.Println(strings.Repeat("-", 80))
		fmt.Printf("%-25s%11s%15s%8s%11s%8s\n", "Sample", "Total", "Non-Mutant", "%", "Mutant", "%")
		fmt.Println(strings.Repeat("-", 80))
		for sample, total := range totalMap {
			fmt.Printf("%-25s%11d%15d%8.2f%11d%8.2f\n", sample, total, countMap[sample], percent(countMap[sample], total), total-countMap[sample], percent(total-countMap[sample], total))
		}

		fmt.Println()
	}
}
