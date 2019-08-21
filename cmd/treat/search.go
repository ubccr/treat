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
	"encoding/csv"
	"fmt"
	"os"
	"strconv"

	"github.com/sirupsen/logrus"
	"github.com/ubccr/treat"
)

func Search(dbpath string, fields *SearchFields, csvOutput, noHeader, fastaOutput bool) {
	s, err := NewStorage(dbpath)
	if err != nil {
		logrus.Fatal(err)
	}

	csvout := csv.NewWriter(os.Stdout)

	if !csvOutput {
		csvout.Comma = '\t'
	}

	if !noHeader && !fastaOutput {
		csvout.Write([]string{
			"gene",
			"sample",
			"norm",
			"read_count",
			"alt_editing",
			"has_mutation",
			"edit_stop",
			"junc_end",
			"junc_len",
			"junc_seq"})
	}

	err = s.Search(fields, func(key *treat.AlignmentKey, a *treat.Alignment) {
		alt := fmt.Sprintf("%d", a.AltEditing)
		if a.AltEditing != 0 {
			alt = fmt.Sprintf("A%d", a.AltEditing)
		}

		if fastaOutput {
			frag, err := s.GetFragment(key, a.Id)
			if err != nil || frag == nil {
				logrus.Printf("fragment not found: %s", err)
				return
			}
			csvout.Write([]string{">" + key.Gene + "|" + key.Sample + "|" + strconv.FormatUint(a.Id, 10) + "|" + frag.Name})
			csvout.Write([]string{frag.String()})
		} else {
			csvout.Write([]string{
				key.Gene,
				key.Sample,
				fmt.Sprintf("%.4f", RoundPlus(a.Norm, 4)),
				fmt.Sprintf("%d", a.ReadCount),
				alt,
				fmt.Sprintf("%d", a.HasMutation),
				fmt.Sprintf("%d", a.EditStop),
				fmt.Sprintf("%d", a.JuncEnd),
				fmt.Sprintf("%d", a.JuncLen),
				a.JuncSeq})
		}

		csvout.Flush()
	})

	if err != nil {
		logrus.Fatal(err)
	}
}
