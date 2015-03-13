package treat

import (
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

    Align(c, tmpl)
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

    grna, err := GrnaFromFile("examples/rps12-guide-rna.csv")
    if err != nil {
        t.Errorf("%s", err)
    }

    aln := NewAlignment(c, tmpl, 0, 0, grna)
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
    if aln.GrnaEditString() != "gRNA13;gRNA14;" {
        t.Errorf("Wrong gnra edit stop. %s != %s", aln.GrnaEditString(), "gRNA13;gRNA14;")
    }
    if aln.GrnaJuncString() != "gRNA14;" {
        t.Errorf("Wrong gnra junc. %s != %s", aln.GrnaJuncString(), "gRNA14;")
    }
}
