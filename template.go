package treat

import (
    "github.com/aebruno/gofasta"
    "os"
    "fmt"
)

type Template struct {
    Bases        string
    EditBase     rune
    EditSite     [][]BaseCountType
}

func NewTemplateFromFasta(path string, orientation OrientationType, base rune) (*Template, error) {
    f, err := os.Open(path)
    if err != nil {
        return nil, fmt.Errorf("Invalid FASTA file: ", err)
    }
    defer f.Close()

    t := make([]*Fragment, 0, 2)

    for rec := range gofasta.SimpleParser(f) {
        frag := NewFragment(rec.Id, rec.Seq, orientation, 0, base)
        t = append(t, frag)
    }

    if len(t) < 2 {
        return nil, fmt.Errorf("Must provide at least 2 templates. Full and Pre edited")
    }

    return NewTemplate(t[0], t[1], t[2:])
}

func NewTemplate(full, pre *Fragment, alt []*Fragment) (*Template, error) {
    if full.Bases != pre.Bases || full.EditBase != pre.EditBase {
        return nil, fmt.Errorf("Invalid template sequences. Full and Pre templates must have the same non-edit bases")
    }

    for _,a := range(alt) {
        if full.Bases != a.Bases || full.EditBase != a.EditBase {
            return nil, fmt.Errorf("Invalid alt template sequence. All templates must have the same non-edit bases")
        }
    }

    editSite := make([][]BaseCountType, len(alt)+2)
    editSite[0] = full.EditSite
    editSite[1] = pre.EditSite
    for i,a := range(alt) {
        editSite[i+2] = a.EditSite
    }

    return &Template{Bases: full.Bases, EditBase: full.EditBase, EditSite: editSite}, nil
}

func (tmpl *Template) Size() (int) {
    return len(tmpl.EditSite)
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
