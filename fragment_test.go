// Copyright 2015 TREAT Authors. All rights reserved.
// Use of this source code is governed by a BSD style
// license that can be found in the LICENSE file.

package treat

import (
    "testing"
    "strings"
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
