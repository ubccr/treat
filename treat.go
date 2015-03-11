package treat

import (
    "github.com/aebruno/nwalgo"
    "github.com/aebruno/gofasta"
    "os"
    "fmt"
    "strings"
    "unicode"
    "bytes"
    "math/big"
)

type OrientationType  int8
type ReadCountType    uint32
type BaseCountType    uint8

const FORWARD OrientationType =   1
const REVERSE OrientationType =  -1

type Alignment struct {
    EditStop       int
    JuncStart      int
    JuncEnd        int
    JuncLen        int
    HasMutation    bool
}

type Template struct {
    Bases        string
    EditBase     rune
    EditSite     [][]BaseCountType
}

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

func NewTemplateFromFasta(path string, orientation OrientationType, base rune) (*Template, error) {
    f, err := os.Open(path)
    if err != nil {
        return nil, fmt.Errorf("Invalid FASTA file: ", err)
    }
    defer f.Close()

    t := make([]*Fragment, 0, 2)

    for rec := range gofasta.SimpleParser(f) {
        frag := NewFragment(rec.Id, rec.Seq, 1, 0, base)
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

func writeBase(buf *bytes.Buffer, base rune, count, max BaseCountType) {
    buf.WriteString(strings.Repeat("-", int(max-count)))
    if count > 0 {
        buf.WriteString(strings.Repeat(string(base), int(count)))
    }
}

func NewAlignment(frag *Fragment, template *Template, primer5, primer3 int) (*Alignment) {
    alignment := new(Alignment)

    m := make([]*big.Int, template.Size())
    for i := range(m) {
        m[i] = new(big.Int)
    }

    aln1, aln2, _ := nwalgo.Align(template.Bases, frag.Bases, 1, -1, -1)

    fi := 0
    ti := 0
    for ai := 0; ai < len(aln1); ai++ {
        if aln1[ai] == '-' {
            fi++
            // insertion
            alignment.HasMutation = true
            continue
        }

        count := BaseCountType(0)
        if aln2[ai] != '-' {
            count = frag.EditSite[fi]

            if frag.Bases[fi] != template.Bases[ti] {
                // SNP
                alignment.HasMutation = true
            }
        } else {
            // deletion
            alignment.HasMutation = true
        }

        for i := range(template.EditSite) {
            if template.EditSite[i][ti] == count {
               m[i].SetBit(m[i], ti, 1)
            }
        }

        if aln2[ai] != '-' {
            fi++
        }
        ti++
    }

    for i := range(template.EditSite) {
        if template.EditSite[i][ti] == frag.EditSite[fi] {
           m[i].SetBit(m[i], ti, 1)
        }
    }

    for j := ti; j >= 0; j-- {
        if m[0].Bit(j) == 0 {
            alignment.JuncStart = ti-j
            break
        }
    }
    for j := 0; j <= ti; j++ {
        if m[1].Bit(j) == 0 {
            alignment.JuncEnd = ti-j
            break
        }
    }

    if alignment.JuncStart > 0 {
        alignment.EditStop = alignment.JuncStart-1
    }

    if alignment.JuncEnd > alignment.JuncStart {
        alignment.JuncLen = alignment.JuncEnd - alignment.EditStop
    }

    return alignment
}

func Align(frag *Fragment, template *Template) {
    aln1, aln2, _ := nwalgo.Align(template.Bases, frag.Bases, 1, -1, -1)

    fragCount := template.Size()+1

    buf := make([]bytes.Buffer, fragCount)
    n := len(aln1)

    fi := 0
    ti := 0
    for ai := 0; ai < n; ai++ {
        if aln1[ai] == '-' {
            for i, _ := range template.EditSite {
                writeBase(&buf[i], template.EditBase, 0, frag.EditSite[fi])
                buf[i].WriteString("-")
            }

            writeBase(&buf[fragCount-1], frag.EditBase, frag.EditSite[fi], frag.EditSite[fi])
            buf[fragCount-1].WriteString(string(frag.Bases[fi]))
            fi++
        } else if aln2[ai] == '-' {
            max := template.Max(ti)

            for i, t := range template.EditSite {
                writeBase(&buf[i], template.EditBase, t[ti], max)
                buf[i].WriteString(string(template.Bases[ti]))
            }
            writeBase(&buf[fragCount-1], '-', 0, max)
            buf[fragCount-1].WriteString("-")
            ti++
        } else {
            max := template.Max(ti)
            if frag.EditSite[fi] > max {
                max = frag.EditSite[fi]
            }

            for i, t := range template.EditSite {
                writeBase(&buf[i], template.EditBase, t[ti], max)
                buf[i].WriteString(string(template.Bases[ti]))
            }
            writeBase(&buf[fragCount-1], frag.EditBase, frag.EditSite[fi], max)
            buf[fragCount-1].WriteString(string(frag.Bases[fi]))
            fi++
            ti++
        }
    }

    // Last edit site has only EditBases
    max := template.Max(ti)
    if frag.EditSite[fi] > max {
        max = frag.EditSite[fi]
    }

    for i, t := range template.EditSite {
        writeBase(&buf[i], template.EditBase, t[ti], max)
    }
    writeBase(&buf[fragCount-1], frag.EditBase, frag.EditSite[fi], max)

    // Write out buffers
    for i,_ := range template.EditSite {
        fmt.Println(buf[i].String())
    }
    fmt.Println(buf[fragCount-1].String())
}
