// Copyright 2015 TREAT Authors. All rights reserved.
// Use of this source code is governed by a BSD style
// license that can be found in the LICENSE file.

package treat

import (
    "strings"
    "unicode"
    "bytes"
    "encoding/gob"
)

type Fragment struct {
    Name         string
    ReadCount    uint64
    Norm         float64
    Bases        string
    EditBase     rune
    EditSite     []BaseCountType
}

// From: http://stackoverflow.com/a/10030772
func reverse(s string) string {
    runes := []rune(s)
    for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
        runes[i], runes[j] = runes[j], runes[i]
    }
    return string(runes)
}

func NewFragment(name, seq string, orientation OrientationType, reads uint64, norm float64, base rune) (*Fragment) {
    base = unicode.ToUpper(base)

    // Ensure all sequences are in forward 5' -> 3' orientation
    if orientation == REVERSE {
        seq = reverse(seq)
    }

    seq = strings.ToUpper(seq)
    base3 := strings.Replace(seq, string(base), "", -1)
    n := len(base3) + 1
    editSite := make([]BaseCountType, n)

    baseCount := BaseCountType(0)
    index := 0
    procBases := func (r rune) rune {
        if r != base {
            editSite[index] = baseCount
            baseCount = 0
            index++
            return r
        }
        baseCount++
        return -1
    }
    bases := strings.Map(procBases, seq)
    editSite[index] = baseCount

    return &Fragment{Name: name, ReadCount: reads, Norm: norm, Bases: bases, EditBase: base, EditSite: editSite}
}

func (f *Fragment) String() (string) {
    var buf bytes.Buffer

    for i,b := range f.EditSite {
        buf.WriteString(strings.Repeat(string(f.EditBase), int(b)))
        if i < len(f.Bases) {
            buf.WriteString(string(f.Bases[i]))
        }
    }

    return buf.String()
}

func (f *Fragment) UnmarshalBytes(data []byte) (error) {
    buf := bytes.NewReader(data)
    dec := gob.NewDecoder(buf)
    err := dec.Decode(&f)
    if err != nil {
        return err
    }

    return nil
}

func (f *Fragment) MarshalBytes() ([]byte, error) {
    data := new(bytes.Buffer)
    enc := gob.NewEncoder(data)
    err := enc.Encode(f)
    if err != nil {
        return nil, err
    }

    return data.Bytes(), nil
}
