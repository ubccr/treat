package treat

import (
    "github.com/aebruno/nwalgo"
    "fmt"
    "strings"
    "bytes"
    "math/big"
)

type Alignment struct {
    EditStop       uint64          `json:"edit_stop"`
    JuncStart      uint64          `json:"junc_start"`
    JuncEnd        uint64          `json:"junc_end"`
    JuncLen        uint64          `json:"junc_len"`
    ReadCount      ReadCountType   `json:"read_count"`
    Norm           float64         `json:"norm_count"`
    HasMutation    uint8           `json:"has_mutation"`
    AltEditing     int8            `json:"alt_editing"`
    GrnaEdit       *big.Int        `json:"-"`
    GrnaJunc       *big.Int        `json:"-"`
}

func writeBase(buf *bytes.Buffer, base rune, count, max BaseCountType) {
    buf.WriteString(strings.Repeat("-", int(max-count)))
    if count > 0 {
        buf.WriteString(strings.Repeat(string(base), int(count)))
    }
}

func NewAlignment(frag *Fragment, template *Template, primer5, primer3 int, grna []*Grna) (*Alignment) {
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
            alignment.HasMutation = uint8(1)
            continue
        }

        count := BaseCountType(0)
        if aln2[ai] != '-' {
            count = frag.EditSite[fi]

            if frag.Bases[fi] != template.Bases[ti] {
                // SNP
                alignment.HasMutation = uint8(1)
            }
        } else {
            // deletion
            alignment.HasMutation = uint8(1)
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

    // Compute binary matrix
    for i := range(template.EditSite) {
        if template.EditSite[i][ti] == frag.EditSite[fi] {
           m[i].SetBit(m[i], ti, 1)
        }
    }

    // Compute junction start site
    for j := ti; j >= 0; j-- {
        if j <= (ti-primer3) && m[0].Bit(j) == 0 {
            alignment.JuncStart = uint64(ti-j)
            break
        }
    }

    // Compute alt editing
    // See if junc start matches an alt template
    alignment.AltEditing = -1
    shift := ti-int(alignment.JuncStart)
    alt := -1
    for i,v := range(m[2:]) {
        if v.Bit(shift) == 1 {
            alt = i
            break
        }
    }

    // If we're at start of alt editing
    if alt > -1 && (ti-shift) == template.AltRegion[alt].Start {
        // Shift Edit Stop Site to first site that doesn't match alt template
        for x := shift; x >= 0; x-- {
            if m[alt+2].Bit(x) == 1 {
                continue
            }

            shift = x
            break
        }

        // If we're before the end of alt editing
        if (ti-shift) > template.AltRegion[alt].End {
            // flag which alt tempalte we matched
            alignment.AltEditing = int8(alt)
            // Shift Junc Start to first site that doesn't match FE template
            for j := shift; j >= 0; j-- {
                if j <= (ti-primer3) && m[0].Bit(j) == 0 {
                    alignment.JuncStart = uint64(ti-j)
                    break
                }
            }
        }
    }

    for j := 0; j <= ti; j++ {
        if j >= primer5 && m[1].Bit(j) == 0 {
            alignment.JuncEnd = uint64(ti-j)
            break
        }
    }

    if alignment.JuncStart > 0 {
        alignment.EditStop = alignment.JuncStart-uint64(1)
    }

    if alignment.JuncEnd > alignment.JuncStart {
        alignment.JuncLen = alignment.JuncEnd - alignment.EditStop
    }

    alignment.ReadCount = frag.ReadCount
    alignment.Norm = frag.Norm


    editVect := big.NewInt(int64(0))
    juncVect := big.NewInt(int64(0))

    // check for gRNA coverage over the edit stop site
    for i, g := range(grna) {
        gstart := g.Start-uint64(11)
        gend := g.End
        sidx := ti-int(alignment.EditStop)
        start := template.BaseIndex[sidx]
        end := template.BaseIndex[sidx]
        if ( (gend >= end && gstart <= start) ||
                (end > gend && start < gstart) ||
                (end <= gstart && end > gend) ||
                (start >= gend && start < gstart) ) {
            editVect.SetBit(editVect, i, 1)
        }
    }

     // check for gRNA coverage over the junction region
    for i, g := range(grna) {
        gstart := g.Start
        gend := g.End
        start := template.BaseIndex[ti-int(alignment.JuncStart)]
        end := template.BaseIndex[ti-int(alignment.JuncEnd)]
        if ( (gend >= end && gstart <= start) ||
                (end > gend && start < gstart) ||
                (end <= gstart && end > gend) ||
                (start >= gend && start < gstart) ) {
            juncVect.SetBit(juncVect, i, 1)
        }
    }

    alignment.GrnaEdit = editVect
    alignment.GrnaJunc = juncVect

    return alignment
}

func (a *Alignment) UnmarshalBinary(data []byte) error {
    var editVect,juncVect uint64

    buf := bytes.NewBuffer(data)
    _, err := fmt.Fscanln(buf,
        &a.EditStop,
        &a.JuncStart,
        &a.JuncEnd,
        &a.JuncLen,
        &a.ReadCount,
        &a.Norm,
        &a.HasMutation,
        &a.AltEditing,
        &editVect,
        &juncVect)

    a.GrnaEdit = big.NewInt(int64(editVect))
    a.GrnaJunc = big.NewInt(int64(juncVect))

    return err
}

func (a *Alignment) MarshalBinary() ([]byte, error) {
    var buf bytes.Buffer
    _, err := fmt.Fprintln(&buf,
        a.EditStop,
        a.JuncStart,
        a.JuncEnd,
        a.JuncLen,
        a.ReadCount,
        a.Norm,
        a.HasMutation,
        a.AltEditing,
        a.GrnaEdit.Uint64(),
        a.GrnaJunc.Uint64())

    return buf.Bytes(), err
}

func (a *Alignment) GrnaEditString() (string) {
    s := ""
    for i := 0; i < a.GrnaEdit.BitLen(); i++ {
        if a.GrnaEdit.Bit(i) == 1 {
            s += fmt.Sprintf("gRNA%d;", i+1)
        }
    }
    return s
}

func (a *Alignment) GrnaJuncString() (string) {
    s := ""
    for i := 0; i < a.GrnaJunc.BitLen(); i++ {
        if a.GrnaJunc.Bit(i) == 1 {
            s += fmt.Sprintf("gRNA%d;", i+1)
        }
    }
    return s
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
