// Copyright 2015 TREAT Authors. All rights reserved.
// Use of this source code is governed by a BSD style
// license that can be found in the LICENSE file.

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
        frag := NewFragment("id", s, FORWARD, 1, 1, 't')

        if strings.ToUpper(s) != frag.String() {
            t.Errorf("%s != %s", s, frag.String())
        }

        // reverse orientation
        frag = NewFragment("id", s, REVERSE, 1, 1, 't')

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

    full := NewFragment("full", "ttCCAATTGCAATTT", FORWARD, 0, 0, 't')
    pre := NewFragment("pre", "ttCAATT", FORWARD, 0, 0, 't')

    _, err = NewTemplate(full, pre, nil, nil)
    if err == nil {
        t.Errorf("Pre and Full templates do not much. Should throw and error")
    }

    full = NewFragment("full", "ttCCAATTGCAATTT", FORWARD, 0, 0, 't')
    pre = NewFragment("pre", "ttttCCAATTTTGCAATTTTT", FORWARD, 0, 0, 't')

    tmpl, err = NewTemplate(full, pre, nil, nil)
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

    fakeAlt := make([]*AltRegion, 2)
    tmpl, err = NewTemplate(full, pre, nil, fakeAlt)
    if err == nil {
        t.Errorf("Alt region should match alt template length. Should throw and error")
    }
}
