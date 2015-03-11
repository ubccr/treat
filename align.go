package treat

import (
    "github.com/aebruno/nwalgo"
    "fmt"
    "strings"
    "bytes"
    "math/big"
    "encoding/binary"
)

type Alignment struct {
    EditStop       uint64
    JuncStart      uint64
    JuncEnd        uint64
    JuncLen        uint64
    ReadCount      ReadCountType
    HasMutation    uint8
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

    for i := range(template.EditSite) {
        if template.EditSite[i][ti] == frag.EditSite[fi] {
           m[i].SetBit(m[i], ti, 1)
        }
    }

    for j := ti; j >= 0; j-- {
        if j <= (ti-primer3) && m[0].Bit(j) == 0 {
            alignment.JuncStart = uint64(ti-j)
            break
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

    return alignment
}

func NewAlignmentFromBytes(data []byte) (*Alignment, error) {
    var a Alignment
    buf := bytes.NewReader(data)

    err := binary.Read(buf, binary.BigEndian, &a.EditStop)
    if err != nil {
        return nil, err
    }
    err = binary.Read(buf, binary.BigEndian, &a.JuncStart)
    if err != nil {
        return nil, err
    }
    err = binary.Read(buf, binary.BigEndian, &a.JuncEnd)
    if err != nil {
        return nil, err
    }
    err = binary.Read(buf, binary.BigEndian, &a.JuncLen)
    if err != nil {
        return nil, err
    }
    err = binary.Read(buf, binary.BigEndian, &a.ReadCount)
    if err != nil {
        return nil, err
    }
    err = binary.Read(buf, binary.BigEndian, &a.HasMutation)
    if err != nil {
        return nil, err
    }

    return &a, nil
}

func (a *Alignment) Bytes() ([]byte, error) {
    data := new(bytes.Buffer)

    err := binary.Write(data, binary.BigEndian, a.EditStop)
    if err != nil {
        return nil, err
    }
    err = binary.Write(data, binary.BigEndian, a.JuncStart)
    if err != nil {
        return nil, err
    }
    err = binary.Write(data, binary.BigEndian, a.JuncEnd)
    if err != nil {
        return nil, err
    }
    err = binary.Write(data, binary.BigEndian, a.JuncLen)
    if err != nil {
        return nil, err
    }
    err = binary.Write(data, binary.BigEndian, a.ReadCount)
    if err != nil {
        return nil, err
    }
    err = binary.Write(data, binary.BigEndian, a.HasMutation)
    if err != nil {
        return nil, err
    }

    return data.Bytes(), nil
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
