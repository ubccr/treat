package treat

import (
    "testing"
    "strings"
    "bytes"
)

func TestFragment(t *testing.T) {
    seqs := []string{"CTGGTTTCA", "CTAATACACTTTTGAttttt", "TtTCTGGCTATT", "ttCGAGTATTT"}
    for _, s := range seqs {
        // forward orientation
        frag := NewFragment("id", s, FORWARD, 1, 't')

        if strings.ToUpper(s) != frag.String() {
            t.Errorf("%s != %s", s, frag.String())
        }

        // reverse orientation
        frag = NewFragment("id", s, REVERSE, 1, 't')

        if strings.ToUpper(reverse(s)) != frag.String() {
            t.Errorf("%s != %s", s, frag.String())
        }
    }
}

func TestTemplate(t *testing.T) {
    tmpl, err := NewTemplateFromFasta("examples/templates.fa", FORWARD, 't')
    if err != nil {
        t.Errorf("%s", err)
    }

    if tmpl.Size() != 5 {
        t.Errorf("Wrong template size. %d != %d", tmpl.Size(), 5)
    }

    full := NewFragment("full", "ttCCAATTGCAATTT", FORWARD, 0, 't')
    pre := NewFragment("pre", "ttCAATT", FORWARD, 0, 't')

    _, err = NewTemplate(full, pre, nil)
    if err == nil {
        t.Errorf("Pre and Full templates do not much. Should throw and error")
    }

    full = NewFragment("full", "ttCCAATTGCAATTT", FORWARD, 0, 't')
    pre = NewFragment("pre", "ttttCCAATTTTGCAATTTTT", FORWARD, 0, 't')

    tmpl, err = NewTemplate(full, pre, nil)
    if err != nil {
        t.Errorf("%s", err)
    }

    var buf bytes.Buffer

    for i,b := range tmpl.EditSite[1] {
        buf.WriteString(strings.Repeat(string(tmpl.EditBase), int(b)))
        if i < len(tmpl.Bases) {
            buf.WriteString(string(tmpl.Bases[i]))
        }
    }

    if pre.String() != buf.String() {
        t.Errorf("%s != %s", full.String(), buf.String())
    }

}

func TestAlign(t *testing.T) {
    fe := "TTTCTGAGTTTAGTAT"
    pe := "TTTTTTCTTTTGAGTTTTTTAGTATT"
    cl := "TTTTCTTTGAGTTTAGTATTT"

    a := NewFragment("a", fe, FORWARD, 1, 't')
    b := NewFragment("b", pe, FORWARD, 1, 't')
    c := NewFragment("c", cl, FORWARD, 1, 't')

    tmpl, err := NewTemplate(a, b, nil)
    if err != nil {
        t.Errorf("%s", err)
    }

    Align(c, tmpl)
}

func TestAlignment(t *testing.T) {
    fe := "CTAATACACTTTTGATAACAAACTAAAGTAAAtAtAttttGttttttttGCGtAtGtGATTTTTGtAtGGttGttGtttACGttttGttttAtttGttttAtGttAttAtAtGAGtCCGCGAttGCCCAGttCCGGtAACCGACGtGtAttGtAtGCCGtAttttAttTAtAtAAttttGtttGGAtGttGCGttGttttttttGttGttttAttGGtttAGttAtGTCAttAtttAttAtAGAGGGTGGtGGttttGttGAtttACCCGGtGTAAAGtAttAtACACGTAttGtAAGttAGATTTAGAtATAAGATATGTTTTT"
    pe := "CTAATACACTTTTGATAACAAACTAAAGTAAAAAGGCGAGGATTTTTTGAGTGGGACTGGAGAGAAAGAGCCGTTCGAGCCCAGCCGGAACCGACGGAGAGCTTCTTTTGAATAAAAGGGAGGCGGGGAGGAGAGTTTCAAAAAGATTTGGGTGGGGGGAACCCTTTGTTTTGGTTAAAGAAACATCGTTTAGAAGAGATTTTAGAATAAGATATGTTTTT"
    cl := "CTAATACACTTTTGATAACAAACTAAAGATATAATATTTTTGTTTTTTTTGCGTATGTGATTTTTGTATGGTTGTTGTTTACGTTTTGTTTTATTTGTTTTATGTTATTATATGAGTCCGCGATTGCCCAGTTCCGGTAACCGACGTGTATTGTATGCCGTATTTTATTTATATAATTTTGTTTGGATGTTGCGTTGTTTTTTTTGTTGTTTTATTGGTTTAGTTATGTCATTATTTATTATAGAGGGTGGTGGTTTTGTTGATTTACCCGGTGTAAAGTATTATACACGTATTGTAAGTTAGATTTAGATATAAGATATGTTTTT"

    a := NewFragment("a", fe, FORWARD, 1, 't')
    b := NewFragment("b", pe, FORWARD, 1, 't')
    c := NewFragment("c", cl, FORWARD, 1, 't')

    tmpl, err := NewTemplate(a, b, nil)
    if err != nil {
        t.Errorf("%s", err)
    }

    aln := NewAlignment(c, tmpl, 0, 0)
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
