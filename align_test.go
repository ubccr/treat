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
	"fmt"
	"os"
	"regexp"
	"strconv"
	"strings"
	"testing"

	"github.com/aebruno/gofasta"
)

func TestAlign(t *testing.T) {
	fe := "TTTCTGAGTTTAGTAT"
	pe := "TTTTTTCTTTTGAGTTTTTTAGTATT"
	cl := "TTTTCTTTGAGTTTAGTATTT"

	a := NewFragment("a", fe, FORWARD, 't')
	b := NewFragment("b", pe, FORWARD, 't')
	c := NewFragment("c", cl, FORWARD, 't')

	tmpl, err := NewTemplate(a, b, nil, nil)
	if err != nil {
		t.Errorf("%s", err)
	}

	buf := new(bytes.Buffer)

	aln := NewAlignment(c, tmpl, false)
	aln.WriteTo(buf, c, tmpl, 80)

	//fmt.Printf("%s", buf.String())
}

func parseKeyVal(id string) map[string]int {
	attr := make(map[string]int)
	pattern := regexp.MustCompile(`[^=\s]+=\-?\d+`)
	matches := pattern.FindAllString(id, -1)

	for _, kv := range matches {
		parts := strings.Split(kv, "=")
		ival, _ := strconv.Atoi(parts[1])
		attr[parts[0]] = ival
	}

	return attr
}

func TestAlignFragments(t *testing.T) {
	tmpl, err := NewTemplateFromFasta("examples/test-templates.fa", FORWARD, 't')
	if err != nil {
		t.Errorf("%s", err)
	}

	if tmpl.Size() != 3 {
		t.Errorf("Wrong template size. %d != %d", tmpl.Size(), 5)
	}

	f, err := os.Open("examples/test-sample.fa")
	if err != nil {
		t.Errorf("Failed to open test sample data")
	}
	defer f.Close()

	for rec := range gofasta.SimpleParser(f) {
		attr := parseKeyVal(rec.Id)
		frag := NewFragment(rec.Id, rec.Seq, FORWARD, rune('t'))
		aln := NewAlignment(frag, tmpl, false)
		if int(aln.EditStop) != attr["ess"] {
			buf := new(bytes.Buffer)
			aln.WriteTo(buf, frag, tmpl, 80)
			fmt.Printf("%s", buf.String())
			t.Errorf("Wrong ESS. %d != %d for sequence id: %s", int(aln.EditStop), attr["ess"], rec.Id)
		}
		if int(aln.JuncStart) != attr["jss"] {
			buf := new(bytes.Buffer)
			aln.WriteTo(buf, frag, tmpl, 80)
			fmt.Printf("%s", buf.String())
			t.Errorf("Wrong JSS. %d != %d for sequence id: %s", int(aln.JuncStart), attr["jss"], rec.Id)
		}
		if int(aln.JuncEnd) != attr["jes"] {
			buf := new(bytes.Buffer)
			aln.WriteTo(buf, frag, tmpl, 80)
			fmt.Printf("%s", buf.String())
			t.Errorf("Wrong JES. %d != %d for sequence id: %s", int(aln.JuncEnd), attr["jes"], rec.Id)
		}
		if int(aln.JuncLen) != (attr["jes"] - attr["ess"]) {
			buf := new(bytes.Buffer)
			aln.WriteTo(buf, frag, tmpl, 80)
			fmt.Printf("%s", buf.String())
			t.Errorf("Wrong Junc Length. %d != %d for sequence id: %s", int(aln.JuncLen), (attr["jes"] - attr["ess"]), rec.Id)
		}
		if int(aln.AltEditing) != attr["alt"] {
			buf := new(bytes.Buffer)
			aln.WriteTo(buf, frag, tmpl, 80)
			fmt.Printf("%s", buf.String())
			t.Errorf("Wrong Alt Editing. %d != %d for sequence id: %s", int(aln.AltEditing), attr["alt"], rec.Id)
		}
		if int(aln.HasMutation) != attr["has_mutation"] {
			buf := new(bytes.Buffer)
			aln.WriteTo(buf, frag, tmpl, 80)
			fmt.Printf("%s", buf.String())
			t.Errorf("Wrong Has Mutation. %d != %d for sequence id: %s", int(aln.HasMutation), attr["has_mutation"], rec.Id)
		}
		if int(aln.Mismatches) != attr["mismatches"] {
			buf := new(bytes.Buffer)
			aln.WriteTo(buf, frag, tmpl, 80)
			fmt.Printf("%s", buf.String())
			t.Errorf("Wrong Mismatch count. %d != %d for sequence id: %s", int(aln.Mismatches), attr["mismatches"], rec.Id)
		}
		if int(aln.Indel) != attr["indel"] {
			buf := new(bytes.Buffer)
			aln.WriteTo(buf, frag, tmpl, 80)
			fmt.Printf("%s", buf.String())
			t.Errorf("Wrong Indel. %d != %d for sequence id: %s", int(aln.Indel), attr["indel"], rec.Id)
		}
	}
}

func BenchmarkBinaryMarshal(b *testing.B) {
	a := new(Alignment)
	x := new(Alignment)
	for i := 0; i < b.N; i++ {
		buf, _ := a.MarshalBinary()
		x.UnmarshalBinary(buf)
	}
}
