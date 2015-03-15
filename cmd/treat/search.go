package main

import (
    "log"
    "fmt"
    "strings"
    "github.com/ubccr/treat"
)


func Search(dbpath string, fields *treat.SearchFields) {
    s, err := treat.NewStorage(dbpath)
    if err != nil {
        log.Fatalf("%s", err)
    }

    err = s.Search(fields, func (k string, a *treat.Alignment) {
        key := strings.Split(string(k), ";")

        alt := fmt.Sprintf("%d", a.AltEditing)
        if a.AltEditing != -1 {
            alt = fmt.Sprintf("A%d", a.AltEditing)
        }

        fmt.Println(strings.Join([]string{
            key[0],
            key[1],
            key[2],
            key[3],
            fmt.Sprintf("%.4f", treat.RoundPlus(a.Norm, 4)),
            fmt.Sprintf("%d", a.ReadCount),
            alt,
            fmt.Sprintf("%d", a.HasMutation),
            fmt.Sprintf("%d", a.EditStop),
            fmt.Sprintf("%d", a.JuncEnd),
            fmt.Sprintf("%d", a.JuncLen),
            a.GrnaEditString()}, "\t"))
    })

    if err != nil {
        log.Fatalf("%s", err)
    }
}
