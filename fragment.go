package treat

import (
    "strings"
    "unicode"
    "bytes"
)

type Fragment struct {
    Name         string
    ReadCount    ReadCountType
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

func NewFragment(name, seq string, orientation OrientationType, reads ReadCountType, base rune) (*Fragment) {
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

    return &Fragment{Name: name, ReadCount: reads, Bases: bases, EditBase: base, EditSite: editSite}
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
