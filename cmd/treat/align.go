package main

import (
    "os"
    "log"
    "bufio"
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

func Align(options *AlignOptions) {
    if len(options.TemplatePath) == 0 {
        log.Fatalln("ERROR Please provide path to templates file")
    }
    if len(options.FragmentPath) == 0 {
        log.Fatalln("ERROR Please provide path to fragment file")
    }
    if len(options.EditBase) != 1 {
        log.Fatalln("ERROR Please provide the edit base")
    }

    grna := make([]*treat.Grna, 0)
    if len(options.GrnaPath) != 0 {
        g, err := treat.GrnaFromFile(options.GrnaPath)
        if err != nil {
            log.Fatalf("ERROR parsing grna file: %s", err)
        }
        grna = g
    }

    tmpl, err := treat.NewTemplateFromFasta(options.TemplatePath, treat.FORWARD, rune(options.EditBase[0]))
    if err != nil {
        log.Fatal(err)
    }

    tmpl.Grna = grna
    tmpl.Primer3 = options.Primer3
    tmpl.Primer5 = options.Primer5

    f, err := os.Open(options.FragmentPath)
    if err != nil {
        log.Fatal(err)
    }
    defer f.Close()

    for rec := range gofasta.SimpleParser(f) {
        frag := treat.NewFragment(rec.Id, rec.Seq, treat.FORWARD, 0, 0, rune(options.EditBase[0]))
        aln := treat.NewAlignment(frag, tmpl)
        buf := bufio.NewWriter(os.Stdout)
        defer buf.Flush()
        aln.WriteTo(buf, frag, tmpl, 80)
    }
}
