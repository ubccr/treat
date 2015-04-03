// Copyright 2015 TREAT Authors. All rights reserved.
// Use of this source code is governed by a BSD style
// license that can be found in the LICENSE file.

package treat

import (
    "github.com/aebruno/gofasta"
    "os"
    "fmt"
    "regexp"
    "strconv"
    "bytes"
    "encoding/gob"
)

var startPattern = regexp.MustCompile(`\s*alt_start=(\d+)\s*`)
var endPattern = regexp.MustCompile(`\s*alt_stop=(\d+)\s*`)

type AltRegion struct {
    Start     int
    End       int
}

type Template struct {
    Bases        string
    EditStop     uint64
    EditBase     rune
    EditSite     [][]BaseCountType
    BaseIndex    []uint64
    AltRegion    []*AltRegion
    Grna         []*Grna
    Primer5      int
    Primer3      int
}

func NewTemplateFromFasta(path string, orientation OrientationType, base rune) (*Template, error) {
    f, err := os.Open(path)
    if err != nil {
        return nil, fmt.Errorf("Invalid FASTA file: %s", err)
    }
    defer f.Close()

    t := make([]*Fragment, 0, 2)
    alt := make([]*AltRegion, 0)

    for rec := range gofasta.SimpleParser(f) {
        frag := NewFragment(rec.Id, rec.Seq, orientation, 0, 0, base)
        t = append(t, frag)

        if len(t) > 2 {
            start, end := -1, -1

            matches := startPattern.FindStringSubmatch(rec.Id)
            if len(matches) == 2 {
                start, err = strconv.Atoi(matches[1])
                if err != nil {
                    start = 0
                }
            }
            matches = endPattern.FindStringSubmatch(rec.Id)
            if len(matches) == 2 {
                end, err = strconv.Atoi(matches[1])
                if err != nil {
                    end = 0
                }
            }

            if start > -1 && end > -1 {
                alt = append(alt, &AltRegion{Start: start, End: end})
            }
        }
    }

    if len(t) < 2 {
        return nil, fmt.Errorf("Must provide at least 2 templates. Full and Pre edited")
    }


    return NewTemplate(t[0], t[1], t[2:], alt)
}

func NewTemplate(full, pre *Fragment, alt []*Fragment, altRegion []*AltRegion) (*Template, error) {
    if full.Bases != pre.Bases || full.EditBase != pre.EditBase {
        return nil, fmt.Errorf("Invalid template sequences. Full and Pre templates must have the same non-edit bases")
    }

    for _,a := range(alt) {
        if full.Bases != a.Bases || full.EditBase != a.EditBase {
            return nil, fmt.Errorf("Invalid alt template sequence. All templates must have the same non-edit bases")
        }
    }

    if len(alt) != len(altRegion) {
        return nil, fmt.Errorf("Invalid alt templates. Please specify the alt regions")
    }

    editSite := make([][]BaseCountType, len(alt)+2)
    editSite[0] = full.EditSite
    editSite[1] = pre.EditSite
    for i,a := range(alt) {
        editSite[i+2] = a.EditSite
    }

    bi := make([]uint64, len(editSite[0]))

    index := uint64(0)
    for i := range(editSite[0]) {
        max := editSite[0][i]
        if editSite[1][i] > max {
            max = editSite[1][i]
        }
        index += uint64(max)
        if i > 0 && i != len(editSite[0])-1 {
            index++
        }
        bi[i] = index
    }

    return &Template{Bases: full.Bases, EditBase: full.EditBase, EditSite: editSite, AltRegion: altRegion, BaseIndex: bi}, nil
}

func (tmpl *Template) Size() (int) {
    return len(tmpl.EditSite)
}

func (tmpl *Template) Len() (int) {
    return len(tmpl.EditSite[0])
}

func (tmpl *Template) Max(i int) (BaseCountType) {
    max := BaseCountType(0)
    for j := range(tmpl.EditSite) {
        if tmpl.EditSite[j][i] > max {
            max = tmpl.EditSite[j][i]
        }
    }
    return max
}

func (tmpl *Template) UnmarshalBytes(data []byte) (error) {
    buf := bytes.NewReader(data)
    dec := gob.NewDecoder(buf)
    err := dec.Decode(&tmpl)
    if err != nil {
        return err
    }

    return nil
}

func (tmpl *Template) MarshalBytes() ([]byte, error) {
    data := new(bytes.Buffer)
    enc := gob.NewEncoder(data)
    err := enc.Encode(tmpl)
    if err != nil {
        return nil, err
    }

    return data.Bytes(), nil
}
