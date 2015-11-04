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

const (
	COUNT_UNIQUE = iota
	COUNT_FRAG
	COUNT_NORM
)

type Stats struct {
	Total          int
	Std            int
	NonStd         int
	Indels         int
	Snps           int
	SingleMismatch int
	DoubleMismatch int
}

type SampleStats struct {
	Stats
}

type GeneStats struct {
	Stats
	Name      string
	SampleMap map[string]*SampleStats
}

func percent(x, y int) float64 {
	if y == 0 {
		return float64(0)
	}

	return (float64(x) / float64(y)) * float64(100)
}

func ShowStats(dbpath, gene string, unique bool, norm bool) {
	s, err := NewStorage(dbpath)
	if err != nil {
		logrus.Fatal(err)
	}

	countby := COUNT_FRAG
	if unique {
		countby = COUNT_UNIQUE
	} else if norm {
		countby = COUNT_NORM
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

		stats, err := geneStats(s, g, countby)
		if err != nil {
			logrus.Fatal(err)
		}

		fmt.Println(strings.Repeat("=", 80))
		fmt.Println(stats.Name)
		fmt.Println(strings.Repeat("=", 80))
		fmt.Printf("%20s%11d\n", "Total Alignments:", stats.Total)
		fmt.Printf("%20s%11d\n", "Standard:", stats.Std)
		fmt.Printf("%20s%11d\n", "Non-Standard:", stats.NonStd)
		fmt.Printf("%20s%11d\n", "1-Mismatch:", stats.SingleMismatch)
		fmt.Printf("%20s%11d\n", "2-Mismatch:", stats.DoubleMismatch)
		fmt.Printf("%20s%11d\n", ">3-Mismatch:", stats.Snps)
		fmt.Printf("%20s%11d\n", "Indels:", stats.Indels)
		fmt.Printf("%20s%11d\n", "Template Edit Stop:", tmpl.EditStop)
		fmt.Printf("%20s%11s\n", "Edit Base:", string(tmpl.EditBase))
		fmt.Printf("%20s%11d\n", "Alt Templates:", len(tmpl.AltRegion))
		fmt.Println(strings.Repeat("-", 80))
		fmt.Printf("%-15s%9s%9s%5s%9s%5s%9s%5s%9s%5s\n", "Sample", "Total", "Std", "%", "Non-Std", "%", "1MM", "%", "2MM", "%")
		fmt.Println(strings.Repeat("-", 80))
		for sample, rec := range stats.SampleMap {
			name := sample
			if len(sample) > 12 {
				name = name[0:12] + ".."
			}
			fmt.Printf("%-15s%9d%9d%5.1f%9d%5.1f%9d%5.1f%9d%5.1f\n",
				name,
				rec.Total,
				rec.Std,
				percent(rec.Std, rec.Total),
				rec.NonStd,
				percent(rec.NonStd, rec.Total),
				rec.SingleMismatch,
				percent(rec.SingleMismatch, rec.Total),
				rec.DoubleMismatch,
				percent(rec.DoubleMismatch, rec.Total))
		}

		fmt.Println()
	}
}

func geneStats(s *Storage, gene string, countby int) (*GeneStats, error) {
	gstat := &GeneStats{Name: gene}
	gstat.SampleMap = make(map[string]*SampleStats)

	err := s.Search(&SearchFields{Gene: gene, All: true, EditStop: -1, JuncLen: -1, JuncEnd: -1}, func(key *treat.AlignmentKey, a *treat.Alignment) {
		if _, ok := gstat.SampleMap[key.Sample]; !ok {
			gstat.SampleMap[key.Sample] = &SampleStats{}
		}

		var readCount int
		if countby == COUNT_NORM {
			readCount = int(a.Norm)
		} else if countby == COUNT_FRAG {
			readCount = int(a.ReadCount)
		} else {
			readCount = 1
		}

		if a.HasMutation == uint8(0) {
			gstat.SampleMap[key.Sample].Std += readCount
			gstat.Std += readCount
		} else {
			gstat.SampleMap[key.Sample].NonStd += readCount
			gstat.NonStd += readCount
		}

		if a.Indel == uint8(1) {
			gstat.SampleMap[key.Sample].Indels += readCount
			gstat.Indels += readCount
		} else if a.Mismatches == uint8(1) {
			gstat.SampleMap[key.Sample].SingleMismatch += readCount
			gstat.SingleMismatch += readCount
		} else if a.Mismatches == uint8(2) {
			gstat.SampleMap[key.Sample].DoubleMismatch += readCount
			gstat.DoubleMismatch += readCount
		} else if a.Mismatches > uint8(2) {
			gstat.SampleMap[key.Sample].Snps += readCount
			gstat.Snps += readCount
		}

		gstat.SampleMap[key.Sample].Total += readCount
		gstat.Total += readCount
	})

	if err != nil {
		return nil, err
	}

	return gstat, nil
}
