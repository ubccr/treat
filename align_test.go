package treat

import (
    "testing"
    "os"
    "math/big"
    "github.com/aebruno/gofasta"
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

    tmpl.Grna = grna

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
    if aln.GrnaEditString() != "gRNA13;gRNA14;" {
        t.Errorf("Wrong grna edit stop. %s != %s", aln.GrnaEditString(), "gRNA13;gRNA14;")
    }
    if aln.GrnaJuncString() != "gRNA14;" {
        t.Errorf("Wrong grna junc. %s != %s", aln.GrnaJuncString(), "gRNA14;")
    }
}

func TestAlignGrna(t *testing.T) {
    grna, err := GrnaFromFile("examples/rps12-guide-rna.csv")
    if err != nil {
        t.Errorf("%s", err)
    }

    tmpl, err := NewTemplateFromFasta("examples/templates.fa", FORWARD, 't')
    if err != nil {
        t.Errorf("%s", err)
    }

    tmpl.Grna = grna

    f, err := os.Open("examples/clones.fa")
    if err != nil {
        t.Errorf("%s", err)
    }
    defer f.Close()

    count := 0
    for rec := range gofasta.SimpleParser(f) {
        frag := NewFragment(rec.Id, rec.Seq, FORWARD, 0, 0, 't')
        aln := NewAlignment(frag, tmpl, 10, 10)
        if count == 0 || count == 1 {
            if aln.GrnaEditString() != "gRNA13;gRNA14;" {
                t.Errorf("Wrong edit grna. %s != %s", aln.GrnaEditString(), "gRNA13;gRNA14;")
            }
        } else if count == 2 {
            if aln.GrnaEditString() != "gRNA7;" {
                t.Errorf("Wrong edit grna. %s != %s", aln.GrnaEditString(), "gRNA7;")
            }
        } else if count == 3 {
            if aln.GrnaEditString() != "gRNA8;gRNA9;" {
                t.Errorf("Wrong edit grna. %s != %s", aln.GrnaEditString(), "gRNA8;gRNA9;")
            }
        } else if count == 4 {
            if aln.GrnaEditString() != "gRNA5;" {
                t.Errorf("Wrong edit grna. %s != %s", aln.GrnaEditString(), "gRNA5;")
            }
        }
        count++
    }
}

func BenchmarkBinaryMarshal(b *testing.B) {
    a := new(Alignment)
    a.GrnaEdit = big.NewInt(int64(0))
    a.GrnaJunc = big.NewInt(int64(0))
    x := new(Alignment)
    for i := 0; i < b.N; i++ {
        buf, _ := a.MarshalBinary()
        x.UnmarshalBinary(buf)
    }
}
