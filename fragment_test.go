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
		frag := NewFragment("id", s, FORWARD, 1, 1, 't')

		if strings.ToUpper(s) != frag.String() {
			t.Errorf("%s != %s", s, frag.String())
		}

		// reverse orientation
		frag = NewFragment("id", s, REVERSE, 1, 1, 't')

		if strings.ToUpper(reverse(s)) != frag.String() {
			t.Errorf("%s != %s", s, frag.String())
		}
	}
}
