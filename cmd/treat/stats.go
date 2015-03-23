package main

import (
    "log"
    "fmt"
    "strings"
    "github.com/ubccr/treat"
)

func percent(x, y int) (float64) {
    return (float64(x)/float64(y))*float64(100)
}

func Stats(dbpath, gene string) {
    s, err := NewStorage(dbpath)
    if err != nil {
        log.Fatalf("%s", err)
    }

    geneTemplates, err = s.TemplateMap()
    if err != nil {
        log.Fatal(err)
    }

    for g, tmpl := range(geneTemplates) {
        if len(gene) > 0 && g != gene {
            continue
        }

        grandTotal := 0
        countMap := make(map[string]int)
        totalMap := make(map[string]int)
        err = s.Search(&SearchFields{Gene: g, All: true, EditStop: -1, JuncLen: -1, JuncEnd: -1}, func (key *treat.AlignmentKey, a *treat.Alignment) {
            if a.HasMutation == uint64(0) {
                countMap[key.Sample]++
            }

            totalMap[key.Sample]++
            grandTotal++
        })

        if err != nil {
            log.Fatalf("%s", err)
        }

        fmt.Println(strings.Repeat("=", 80))
        fmt.Println(g)
        fmt.Println(strings.Repeat("=", 80))
        fmt.Printf("%20s%11d\n", "Total Alignments:", grandTotal)
        fmt.Printf("%20s%11d\n", "Alt Templates:", len(tmpl.Grna))
        fmt.Printf("%20s%11d\n", "Guide RNAs:", len(tmpl.AltRegion))
        fmt.Println(strings.Repeat("-", 80))
        fmt.Printf("%-30s%11s%11s%8s%11s%8s\n", "Sample", "Total", "OK", "%", "Mutant", "%")
        fmt.Println(strings.Repeat("-", 80))
        for sample, total := range(totalMap) {
            fmt.Printf("%-30s%11d%11d%8.2f%11d%8.2f\n", sample, total, countMap[sample], percent(countMap[sample], total), total-countMap[sample], percent(total-countMap[sample], total))
        }

        fmt.Println()
    }
}
