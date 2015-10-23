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
	"bytes"
	"strings"
	"testing"
)

func TestTemplate(t *testing.T) {
	tmpl, err := NewTemplateFromFasta("examples/templates.fa", FORWARD, 't')
	if err != nil {
		t.Errorf("%s", err)
	}

	if tmpl.Size() != 5 {
		t.Errorf("Wrong template size. %d != %d", tmpl.Size(), 5)
	}

	err = tmpl.SetPrimer5("XXX")
	if err == nil {
		t.Errorf("Setting invald primer region should throw error")
	}

	err = tmpl.SetPrimer5("CTAATACACTTTTGATAACAAAC")
	if err != nil {
		t.Errorf("%s", err)
	}

	if tmpl.Primer5 != 16 {
		t.Errorf("Primer5 region should be 16 edit sites got: %d", tmpl.Primer5)
	}

	err = tmpl.SetPrimer3("XXX")
	if err == nil {
		t.Errorf("Setting invald primer region should throw error")
	}

	err = tmpl.SetPrimer3("AAAAACATATCTTATATCTAAA")
	if err != nil {
		t.Errorf("%s", err)
	}

	if tmpl.Primer3 != 10 {
		t.Errorf("Primer3 region should be 10 edit sites got: %d", tmpl.Primer3)
	}

	full := NewFragment("full", "ttCCAATTGCAATTT", FORWARD, 0, 0, 't')
	pre := NewFragment("pre", "ttCAATT", FORWARD, 0, 0, 't')

	_, err = NewTemplate(full, pre, nil, nil)
	if err == nil {
		t.Errorf("Pre and Full templates do not much. Should throw and error")
	}

	full = NewFragment("full", "ttCCAATTGCAATTT", FORWARD, 0, 0, 't')
	pre = NewFragment("pre", "ttttCCAATTTTGCAATTTTT", FORWARD, 0, 0, 't')

	tmpl, err = NewTemplate(full, pre, nil, nil)
	if err != nil {
		t.Errorf("%s", err)
	}

	var buf bytes.Buffer

	for i, b := range tmpl.EditSite[1] {
		buf.WriteString(strings.Repeat(string(tmpl.EditBase), int(b)))
		if i < len(tmpl.Bases) {
			buf.WriteString(string(tmpl.Bases[i]))
		}
	}

	if pre.String() != buf.String() {
		t.Errorf("%s != %s", full.String(), buf.String())
	}

	fakeAlt := make([]*AltRegion, 2)
	tmpl, err = NewTemplate(full, pre, nil, fakeAlt)
	if err == nil {
		t.Errorf("Alt region should match alt template length. Should throw and error")
	}
}
