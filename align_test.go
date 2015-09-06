// Copyright 2015 TREAT Authors. All rights reserved.
// Use of this source code is governed by a BSD style
// license that can be found in the LICENSE file.

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
