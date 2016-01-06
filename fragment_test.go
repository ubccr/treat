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

package treat

import (
	"strings"
	"testing"
)

func TestFragment(t *testing.T) {
	seqs := []string{"CTGGTTTCA", "CTAATACACTTTTGAttttt", "TtTCTGGCTATT", "ttCGAGTATTT"}
	for _, s := range seqs {
		// forward orientation
		frag := NewFragment("id", s, FORWARD, 't')

		if strings.ToUpper(s) != frag.String() {
			t.Errorf("%s != %s", s, frag.String())
		}

		if uint32(1) != frag.ReadCount {
			t.Errorf("%d != %d", 1, frag.ReadCount)
		}

		// reverse orientation
		frag = NewFragment("1-143", s, REVERSE, 't')

		if strings.ToUpper(reverse(s)) != frag.String() {
			t.Errorf("%s != %s", s, frag.String())
		}

		if uint32(143) != frag.ReadCount {
			t.Errorf("%d != %d", 143, frag.ReadCount)
		}
	}
}

func TestParseReadCounts(t *testing.T) {
	seqs := map[string]int{
		" 132-2082":                        2082,
		"SAMPLE1_GENE_123432_2082":         2082,
		"132-2082":                         2082,
		"GENE_88772-2082":                  2082,
		"791-2082":                         2082,
		"791-2082            ":             2082,
		"791-2082 attr1=x attr2=y attr3=z": 2082,
		"7912082 attr1=x attr2=y attr3=z":  1,
		"badtest2082":                      1,
		"2082":                             1,
		"-2082":                            1,
		"_2082":                            1,
		"2082_":                            1,
	}

	for id, count := range seqs {
		frag := NewFragment(id, "CTGCTG", FORWARD, 't')

		if uint32(count) != frag.ReadCount {
			t.Errorf("ID %s : %d != %d", id, count, frag.ReadCount)
		}
	}
}
