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
	"testing"
)

func TestAlign(t *testing.T) {
	fe := "TTTCTGAGTTTAGTAT"
	pe := "TTTTTTCTTTTGAGTTTTTTAGTATT"
	cl := "TTTTCTTTGAGTTTAGTATTT"

	a := NewFragment("a", fe, FORWARD, 1, 1, 't')
	b := NewFragment("b", pe, FORWARD, 1, 1, 't')
	c := NewFragment("c", cl, FORWARD, 1, 1, 't')

	tmpl, err := NewTemplate(a, b, nil, nil)
	if err != nil {
		t.Errorf("%s", err)
	}

	buf := new(bytes.Buffer)

	aln := NewAlignment(c, tmpl, false)
	aln.WriteTo(buf, c, tmpl, 80)

	fmt.Printf("%s", buf.String())
}

func TestAlignment(t *testing.T) {
	fe := "CTAATACACTTTTGATAACAAACTAAAGTAAAtAtAttttGttttttttGCGtAtGtGATTTTTGtAtGGttGttGtttACGttttGttttAtttGttttAtGttAttAtAtGAGtCCGCGAttGCCCAGttCCGGtAACCGACGtGtAttGtAtGCCGtAttttAttTAtAtAAttttGtttGGAtGttGCGttGttttttttGttGttttAttGGtttAGttAtGTCAttAtttAttAtAGAGGGTGGtGGttttGttGAtttACCCGGtGTAAAGtAttAtACACGTAttGtAAGttAGATTTAGAtATAAGATATGTTTTT"
	pe := "CTAATACACTTTTGATAACAAACTAAAGTAAAAAGGCGAGGATTTTTTGAGTGGGACTGGAGAGAAAGAGCCGTTCGAGCCCAGCCGGAACCGACGGAGAGCTTCTTTTGAATAAAAGGGAGGCGGGGAGGAGAGTTTCAAAAAGATTTGGGTGGGGGGAACCCTTTGTTTTGGTTAAAGAAACATCGTTTAGAAGAGATTTTAGAATAAGATATGTTTTT"
	cl := "CTAATACACTTTTGATAACAAACTAAAGATATAATATTTTTGTTTTTTTTGCGTATGTGATTTTTGTATGGTTGTTGTTTACGTTTTGTTTTATTTGTTTTATGTTATTATATGAGTCCGCGATTGCCCAGTTCCGGTAACCGACGTGTATTGTATGCCGTATTTTATTTATATAATTTTGTTTGGATGTTGCGTTGTTTTTTTTGTTGTTTTATTGGTTTAGTTATGTCATTATTTATTATAGAGGGTGGTGGTTTTGTTGATTTACCCGGTGTAAAGTATTATACACGTATTGTAAGTTAGATTTAGATATAAGATATGTTTTT"

	a := NewFragment("a", fe, FORWARD, 1, 1, 't')
	b := NewFragment("b", pe, FORWARD, 1, 1, 't')
	c := NewFragment("c", cl, FORWARD, 1, 1, 't')

	tmpl, err := NewTemplate(a, b, nil, nil)
	if err != nil {
		t.Errorf("%s", err)
	}

	aln := NewAlignment(c, tmpl, false)
	if aln.JuncLen != 6 {
		t.Errorf("%s", err)
		t.Errorf("Wrong junc len. %d != %d", aln.JuncLen, 6)
	}
	if aln.EditStop != 137 {
		t.Errorf("%s", err)
		t.Errorf("Wrong edit stop. %d != %d", aln.EditStop, 137)
	}
	if aln.JuncEnd != 143 {
		t.Errorf("%s", err)
		t.Errorf("Wrong edit stop. %d != %d", aln.EditStop, 143)
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
