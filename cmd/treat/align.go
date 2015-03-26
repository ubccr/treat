package main

import (
    "os"
    "log"
    "bufio"
    "fmt"
    "github.com/aebruno/gofasta"
    "github.com/ubccr/treat"
)

type AlignOptions struct {
    TemplatePath  string
    FragmentPath  string
    Primer5       int
    Primer3       int
    GrnaPath      string
    EditBase      string
}

func PrintAlignment(a1, a2 string, tw int) {
    n := len(a1)
    cols := tw-4
    rows := n/cols
    if n % cols > 0 {
        rows++
    }

    for r := 0; r < rows; r++ {
        end := (r*cols)+cols
        if end > n {
            end = n
        }

        fmt.Print("F1: ")
        fmt.Println(a1[(r*cols):end])
        fmt.Print("F2: ")
        fmt.Println(a2[(r*cols):end])
        fmt.Println()
    }
}

func Align(options *AlignOptions) {
    var tmpl *treat.Template

    if len(options.FragmentPath) == 0 {
        log.Fatalln("ERROR Please provide path to fragment file")
    }
    if len(options.EditBase) != 1 {
        log.Fatalln("ERROR Please provide the edit base")
    }

    if len(options.TemplatePath) > 0 {
        t, err := treat.NewTemplateFromFasta(options.TemplatePath, treat.FORWARD, rune(options.EditBase[0]))
        if err != nil {
            log.Fatal(err)
        }
        tmpl = t

        grna := make([]*treat.Grna, 0)
        if len(options.GrnaPath) != 0 {
            g, err := treat.GrnaFromFile(options.GrnaPath)
            if err != nil {
                log.Fatalf("ERROR parsing grna file: %s", err)
            }
            grna = g
        }


        tmpl.Grna = grna
        tmpl.Primer3 = options.Primer3
        tmpl.Primer5 = options.Primer5
    }

    f, err := os.Open(options.FragmentPath)
    if err != nil {
        log.Fatal(err)
    }
    defer f.Close()

    if tmpl == nil {
        frags := make([]*treat.Fragment, 0)
        for rec := range gofasta.SimpleParser(f) {
            frag := treat.NewFragment(rec.Id, rec.Seq, treat.FORWARD, 0, 0, rune(options.EditBase[0]))
            frags = append(frags, frag)
            if len(frags) >= 2 {
                break
            }
        }

        if len(frags) < 2 {
            log.Fatalf("ERROR need 2 fragments to align")
        }

        aln := new(treat.Alignment)
        a1, a2 := aln.SimpleAlign(frags[0], frags[1])
        PrintAlignment(a1, a2, 80)

    } else {
        for rec := range gofasta.SimpleParser(f) {
            frag := treat.NewFragment(rec.Id, rec.Seq, treat.FORWARD, 0, 0, rune(options.EditBase[0]))
            aln := treat.NewAlignment(frag, tmpl)
            buf := bufio.NewWriter(os.Stdout)
            defer buf.Flush()
            aln.WriteTo(buf, frag, tmpl, 80)
        }
    }
}
