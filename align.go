// Copyright 2015 TREAT Authors. All rights reserved.
//
// This file is part of TREAT.
//
// TREAT is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// TREAT is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with TREAT.  If not, see <http://www.gnu.org/licenses/>.

package treat

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"io"
	"math"
	"math/big"
	"strings"

	"github.com/aebruno/nwalgo"
)

type AlignmentKey struct {
	Gene   string
	Sample string
}

type Alignment struct {
	Key         *AlignmentKey `json:"-"`
	Id          uint64        `json:"-"`
	EditStop    uint32        `json:"edit_stop"`
	JuncStart   uint32        `json:"junc_start"`
	JuncEnd     uint32        `json:"junc_end"`
	JuncLen     uint32        `json:"junc_len"`
	ReadCount   uint32        `json:"read_count"`
	Norm        float64       `json:"norm_count"`
	HasMutation uint8         `json:"has_mutation"`
	Mismatches  uint8         `json:"mismatches"`
	Indel       uint8         `json:"indel"`
	AltEditing  uint8         `json:"alt_editing"`
	JuncSeq     string        `json:"-"`
}

func (k *AlignmentKey) UnmarshalBinary(data []byte) error {
	parts := strings.Split(string(data), ";")
	k.Gene = parts[0]
	k.Sample = parts[1]

	return nil
}

func (k *AlignmentKey) MarshalBinary() ([]byte, error) {
	return []byte(strings.Join([]string{k.Gene, k.Sample}, ";")), nil
}

func writeBase(buf *bytes.Buffer, base rune, count, max BaseCountType) {
	buf.WriteString(strings.Repeat("-", int(max-count)))
	if count > 0 {
		buf.WriteString(strings.Repeat(string(base), int(count)))
	}
}

func NewAlignment(frag *Fragment, template *Template, excludeSnps bool) *Alignment {
	alignment := new(Alignment)

	m := make([]*big.Int, template.Size())
	for i := range m {
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
			alignment.Indel = uint8(1)
			continue
		}

		count := BaseCountType(0)
		if aln2[ai] != '-' {
			count = frag.EditSite[fi]

			if frag.Bases[fi] != template.Bases[ti] {
				// SNP
				alignment.Mismatches++

				//TODO: make the number of mismatches configurable
				if excludeSnps || alignment.Mismatches > 2 {
					alignment.HasMutation = uint8(1)
				}
			}
		} else {
			// deletion
			alignment.HasMutation = uint8(1)
			alignment.Indel = uint8(1)
		}

		for i := range template.EditSite {
			if template.EditSite[i][ti] == count {
				m[i].SetBit(m[i], ti, 1)
			}
		}

		if aln2[ai] != '-' {
			fi++
		}
		ti++
	}

	// Last edit site
	for i := range template.EditSite {
		if template.EditSite[i][ti] == frag.EditSite[fi] {
			m[i].SetBit(m[i], ti, 1)
		}
	}

	// Compute junction start site
	isFullyEdited := true
	for j := ti; j >= 0; j-- {
		if m[0].Bit(j) == 0 {
			// Sites are numbered 3' -> 5', so we reverse the index (ti-j)
			alignment.JuncStart = uint32(ti - j)
			isFullyEdited = false
			break
		}
	}

	// If fragment matches the fully edited template entirely.
	if isFullyEdited {
		alignment.JuncStart = uint32(ti)
	}

	// Compute alt editing
	// See if junc start matches an alt template
	shift := ti - int(alignment.JuncStart)
	alt := 0
	hasAlt := false
	for i, v := range m[2:] {
		if v.Bit(shift) == 1 {
			alt = i
			hasAlt = true
			break
		}
	}

	// If we're at start of alt editing
	if hasAlt && (ti-shift) == template.AltRegion[alt].Start {
		// Shift Edit Stop Site to first site that doesn't match alt template
		for x := shift; x >= 0; x-- {
			if m[alt+2].Bit(x) == 1 {
				continue
			}

			shift = x
			break
		}

		// If we're before the end of alt editing
		if (ti - shift) > template.AltRegion[alt].End {
			// flag which alt tempalte we matched
			alignment.AltEditing = uint8(alt + 1)
			// Shift Junc Start to first site that doesn't match FE template
			for j := shift; j >= 0; j-- {
				if j <= ti && m[0].Bit(j) == 0 {
					alignment.JuncStart = uint32(ti - j)
					break
				}
			}
		}
	}

	for j := 0; j <= ti; j++ {
		if m[1].Bit(j) == 0 {
			alignment.JuncEnd = uint32(ti - j)
			break
		}
	}

	if alignment.JuncStart > 0 {
		alignment.EditStop = alignment.JuncStart - uint32(1)
	}

	if alignment.JuncEnd > alignment.EditStop {
		alignment.JuncLen = alignment.JuncEnd - alignment.EditStop
		if alignment.HasMutation == 0 {
			for i := ti - int(alignment.JuncEnd); i < ti-int(alignment.EditStop); i++ {
				alignment.JuncSeq += strings.Repeat(string(frag.EditBase), int(frag.EditSite[i]))
				if i < len(frag.Bases) {
					alignment.JuncSeq += string(frag.Bases[i])
				}
			}
		}
	} else if alignment.JuncEnd < alignment.EditStop {
		alignment.JuncEnd = alignment.EditStop
	}

	// If fragment matches the fully edited template entirely. Set junc len 0
	if isFullyEdited {
		alignment.JuncEnd = alignment.JuncStart
		alignment.JuncLen = 0
		alignment.JuncSeq = ""
	}

	alignment.ReadCount = frag.ReadCount
	alignment.Norm = frag.Norm

	return alignment
}

func (a *Alignment) UnmarshalBinary(buf []byte) error {
	a.EditStop = binary.BigEndian.Uint32(buf[0:4])
	a.JuncStart = binary.BigEndian.Uint32(buf[4:8])
	a.JuncEnd = binary.BigEndian.Uint32(buf[8:12])
	a.JuncLen = binary.BigEndian.Uint32(buf[12:16])
	a.ReadCount = binary.BigEndian.Uint32(buf[16:20])
	a.HasMutation = buf[20]
	a.AltEditing = buf[21]
	a.Mismatches = buf[22]
	a.Indel = buf[23]
	normBits := binary.BigEndian.Uint64(buf[24:32])
	seqLen := binary.BigEndian.Uint32(buf[32:36])

	a.Norm = math.Float64frombits(normBits)

	a.JuncSeq = string(buf[36 : 36+int(seqLen)])

	return nil
}

func (a *Alignment) MarshalBinary() ([]byte, error) {
	buf := make([]byte, 36)

	binary.BigEndian.PutUint32(buf[0:4], a.EditStop)
	binary.BigEndian.PutUint32(buf[4:8], a.JuncStart)
	binary.BigEndian.PutUint32(buf[8:12], a.JuncEnd)
	binary.BigEndian.PutUint32(buf[12:16], a.JuncLen)
	binary.BigEndian.PutUint32(buf[16:20], a.ReadCount)
	buf[20] = a.HasMutation
	buf[21] = a.AltEditing
	buf[22] = a.Mismatches
	buf[23] = a.Indel
	binary.BigEndian.PutUint64(buf[24:32], math.Float64bits(a.Norm))

	seq := []byte(a.JuncSeq)
	binary.BigEndian.PutUint32(buf[32:36], uint32(len(seq)))
	buf = append(buf, seq...)

	return buf, nil
}

func (a *Alignment) SimpleAlign(f1, f2 *Fragment) (string, string) {
	aln1, aln2, _ := nwalgo.Align(f1.Bases, f2.Bases, 1, -1, -1)

	buf := make([]bytes.Buffer, 2)
	n := len(aln1)

	fi := 0
	ti := 0
	for ai := 0; ai < n; ai++ {
		if aln1[ai] == '-' {
			writeBase(&buf[0], f1.EditBase, 0, f2.EditSite[fi])
			buf[0].WriteString("-")

			writeBase(&buf[1], f2.EditBase, f2.EditSite[fi], f2.EditSite[fi])
			buf[1].WriteString(string(f2.Bases[fi]))
			fi++
		} else if aln2[ai] == '-' {
			writeBase(&buf[0], f1.EditBase, f1.EditSite[ti], f1.EditSite[ti])
			buf[0].WriteString(string(f1.Bases[ti]))

			writeBase(&buf[1], '-', 0, f1.EditSite[ti])
			buf[1].WriteString("-")
			ti++
		} else {
			max := f1.EditSite[ti]
			if f2.EditSite[fi] > max {
				max = f2.EditSite[fi]
			}

			writeBase(&buf[0], f1.EditBase, f1.EditSite[ti], max)
			buf[0].WriteString(string(f1.Bases[ti]))

			writeBase(&buf[1], f2.EditBase, f2.EditSite[fi], max)
			buf[1].WriteString(string(f2.Bases[fi]))
			fi++
			ti++
		}
	}

	// Last edit site has only EditBases
	max := f1.EditSite[ti]
	if f2.EditSite[fi] > max {
		max = f2.EditSite[fi]
	}

	writeBase(&buf[0], f1.EditBase, f1.EditSite[ti], max)
	writeBase(&buf[1], f2.EditBase, f2.EditSite[fi], max)

	return buf[0].String(), buf[1].String()
}

func (a *Alignment) WriteTo(w io.Writer, frag *Fragment, template *Template, tw int) error {
	if tw <= 0 {
		tw = 80
	}

	aln1, aln2, _ := nwalgo.Align(template.Bases, frag.Bases, 1, -1, -1)

	fragCount := template.Size() + 1

	buf := make([]bytes.Buffer, fragCount)
	n := len(aln1)

	fi := 0
	ti := 0
	for ai := 0; ai < n; ai++ {
		if aln1[ai] == '-' {
			for i := range template.EditSite {
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

	cols := tw - 4
	rows := buf[0].Len() / cols
	if buf[0].Len()%cols > 0 {
		rows++
	}

	labels := []string{"FE", "PE"}
	for i := range template.AltRegion {
		labels = append(labels, fmt.Sprintf("A%d", i+1))
	}
	labels = append(labels, "CL")

	_, err := w.Write([]byte(strings.Repeat("=", tw) + "\n"))
	if err != nil {
		return err
	}

	_, err = w.Write([]byte(frag.Name + "\n"))
	if err != nil {
		return err
	}

	_, err = w.Write([]byte(strings.Repeat("=", tw) + "\n"))
	if err != nil {
		return err
	}

	_, err = w.Write([]byte(fmt.Sprintf("EditStop: %d\nJuncEnd: %d\nJuncLen: %d\n", a.EditStop, a.JuncEnd, a.JuncLen)))
	if err != nil {
		return err
	}
	_, err = w.Write([]byte(strings.Repeat("=", tw) + "\n\n"))
	if err != nil {
		return err
	}

	for r := 0; r < rows; r++ {
		for i, b := range buf {
			_, err := w.Write([]byte(labels[i] + ": "))
			if err != nil {
				return err
			}

			end := (r * cols) + cols
			if end > b.Len() {
				end = b.Len()
			}
			_, err = w.Write(b.Bytes()[(r * cols):end])
			if err != nil {
				return err
			}
			_, err = w.Write([]byte("\n"))
			if err != nil {
				return err
			}
		}
		_, err := w.Write([]byte("\n"))
		if err != nil {
			return err
		}
	}

	return nil
}
