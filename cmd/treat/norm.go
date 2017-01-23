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
	"github.com/Sirupsen/logrus"
	"github.com/ubccr/treat"
)

func Normalize(dbpath, gene string, norm float64) {
	s, err := NewStorageWrite(dbpath)
	if err != nil {
		logrus.Fatal(err)
	}

	genes, err := s.Genes()
	if err != nil {
		logrus.Fatal(err)
	}

	for _, g := range genes {
		if len(gene) > 0 && g != gene {
			continue
		}
		logrus.Printf("Processing gene %s...", g)

		samples, err := s.SampleKeys(g)
		if err != nil {
			logrus.Fatal(err)
		}

		gnorm := norm

		// Default norm to average read count
		if gnorm == 0 {
			logrus.Info("Using default option of normalizing to average read count across all samples")

			total := 0
			err := s.Search(&SearchFields{Gene: g, HasMutation: false, EditStop: -1, JuncLen: -1, JuncEnd: -1}, func(key *treat.AlignmentKey, a *treat.Alignment) {
				total += int(a.ReadCount)
			})
			if err != nil {
				logrus.Fatal(err)
			}

			gnorm = float64(total) / float64(len(samples))
			logrus.Printf("Total standard reads across all samples: %d", total)
		}

		logrus.Printf("Normalizing to read count: %.4f", gnorm)
		for _, skey := range samples {
			err = s.NormalizeSample(skey, gnorm)
			if err != nil {
				logrus.Fatal(err)
			}
		}
	}
}
