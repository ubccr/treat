// Copyright 2015 TREAT Authors. All rights reserved.
// Use of this source code is governed by a BSD style
// license that can be found in the LICENSE file.

package main

import (
    "log"
    "os"
    "fmt"
    "encoding/csv"
    "github.com/ubccr/treat"
)


func Search(dbpath string, fields *SearchFields, csvOutput, noHeader bool) {
    s, err := NewStorage(dbpath)
    if err != nil {
        log.Fatalf("%s", err)
    }

    csvout := csv.NewWriter(os.Stdout)

    if !csvOutput {
        csvout.Comma = '\t'
    }

    if !noHeader {
        csvout.Write([]string{
            "gene",
            "sample",
            "norm",
            "read_count",
            "alt_editing",
            "has_mutation",
            "edit_stop",
            "junc_end",
            "junc_len",
            "junc_seq"})
    }

    err = s.Search(fields, func (key *treat.AlignmentKey, a *treat.Alignment) {
        alt := fmt.Sprintf("%d", a.AltEditing)
        if a.AltEditing != 0 {
            alt = fmt.Sprintf("A%d", a.AltEditing)
        }

        csvout.Write([]string{
            key.Gene,
            key.Sample,
            fmt.Sprintf("%.4f", RoundPlus(a.Norm, 4)),
            fmt.Sprintf("%d", a.ReadCount),
            alt,
            fmt.Sprintf("%d", a.HasMutation),
            fmt.Sprintf("%d", a.EditStop),
            fmt.Sprintf("%d", a.JuncEnd),
            fmt.Sprintf("%d", a.JuncLen),
            a.JuncSeq})

        csvout.Flush()
    })

    if err != nil {
        log.Fatalf("%s", err)
    }
}
