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
	"strings"

	"github.com/aebruno/nwalgo"
	"github.com/willf/bitset"
)

type AlignmentKey struct {
	Gene   string
	Sample string
}

type Alignment struct {
	Key         *AlignmentKey `json:"-"`
	Id          uint64        `json:"-"`
	EditStop    int           `json:"edit_stop"`
	JuncStart   int           `json:"junc_start"`
	JuncEnd     int           `json:"junc_end"`
	JuncLen     int           `json:"junc_len"`
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

func writeBase(buf *bytes.Buffer, base rune, count, max uint32) {
	buf.WriteString(strings.Repeat("-", int(max-count)))
	if count > 0 {
		buf.WriteString(strings.Repeat(string(base), int(count)))
	}
}

func (a *Alignment) findJES(b *bitset.BitSet) int {
	if b.All() {
		return -1
	} else if b.None() {
		return int(b.Len() - 1)
	} else {
		xor := b.SymmetricDifference(bitset.New(b.Len()).Complement())
		max := uint(0)
		for i, e := xor.NextSet(0); e; i, e = xor.NextSet(i + 1) {
			max = i
		}
		return int(max)
	}
}

func (a *Alignment) findJSS(b *bitset.BitSet) int {
	if b.All() {
		return int(b.Len() - 1)
	} else if b.None() {
		return -1
	} else {
		xor := b.SymmetricDifference(bitset.New(b.Len()).Complement())
		min, _ := xor.NextSet(0)
		return int(min)
	}
}

func (a *Alignment) computeT(frag *Fragment, tmpl *Template, excludeSnps bool) []*bitset.BitSet {
	size := uint(tmpl.Len())
	T := make([]*bitset.BitSet, tmpl.Size())
	for i := range T {
		T[i] = bitset.New(size)
	}

	aln1, aln2, _ := nwalgo.Align(tmpl.Bases, frag.Bases, 1, -1, -1)

	fi := 0
	ti := 0
	for ai := 0; ai < len(aln1); ai++ {
		if aln1[ai] == '-' {
			fi++
			// insertion
			a.HasMutation = uint8(1)
			a.Indel = uint8(1)
			continue
		}

		count := uint32(0)
		if aln2[ai] != '-' {
			count = frag.EditSite[fi]

			if frag.Bases[fi] != tmpl.Bases[ti] {
				// SNP
				a.Mismatches++

				//TODO: make the number of mismatches configurable
				if excludeSnps || a.Mismatches > 2 {
					a.HasMutation = uint8(1)
				}
			}
		} else {
			// deletion
			a.HasMutation = uint8(1)
			a.Indel = uint8(1)
		}

		for i := range tmpl.EditSite {
			if tmpl.EditSite[i][ti] == count {
				T[i] = T[i].Set((size - 1) - uint(ti))
			}
		}

		if aln2[ai] != '-' {
			fi++
		}
		ti++
	}

	// Last edit site
	for i := range tmpl.EditSite {
		if tmpl.EditSite[i][ti] == frag.EditSite[fi] {
			T[i] = T[i].Set((size - 1) - uint(ti))
		}
	}

	return T
}

func (a *Alignment) computeAltEditing(tmpl *Template, T []*bitset.BitSet) {
	// Compute alt editing
	// See if junc start matches an alt template
	shift := a.JuncStart
	alt := 0
	hasAlt := false
	for i, v := range T[2:] {
		if v.Test(uint(shift)) {
			alt = i
			hasAlt = true
			break
		}
	}

	if !hasAlt {
		return
	}

	// If we're not at the start of alt editing return
	if shift != tmpl.AltRegion[alt].Start {
		return

	}

	// Shift Edit Stop Site to first site that doesn't match alt template
	fullMatch := true
	for x := shift; x < tmpl.Len(); x++ {
		if T[alt+2].Test(uint(x)) {
			continue
		}

		fullMatch = false
		shift = x
		break
	}

	if fullMatch || shift >= tmpl.AltRegion[alt].End {
		// If we're after the end of alt editing
		// flag which alt tempalte we matched
		a.AltEditing = uint8(alt + 1)
		// Shift Junc Start to first site that doesn't match FE template
		if fullMatch {
			a.JuncStart = tmpl.Len() - 1
		} else {
			for j := shift; j < tmpl.Len(); j++ {
				if !T[0].Test(uint(j)) {
					a.JuncStart = j
					break
				}
			}
		}
	}
}

func NewAlignment(frag *Fragment, tmpl *Template, excludeSnps bool) *Alignment {
	a := new(Alignment)

	T := a.computeT(frag, tmpl, excludeSnps)

	a.JuncStart = a.findJSS(T[0])
	a.computeAltEditing(tmpl, T)
	a.JuncEnd = a.findJES(T[1])
	a.EditStop = a.JuncStart - 1

	if a.JuncEnd > a.EditStop {
		a.JuncLen = a.JuncEnd - a.EditStop
		if a.HasMutation == 0 {
			from := (tmpl.Len() - 1) - a.JuncEnd
			to := (tmpl.Len() - 1) - a.EditStop
			for i := from; i < to; i++ {
				a.JuncSeq += strings.Repeat(string(frag.EditBase), int(frag.EditSite[i]))
				if i < len(frag.Bases) {
					a.JuncSeq += string(frag.Bases[i])
				}
			}
		}
	} else if a.JuncEnd < a.EditStop {
		a.JuncEnd = a.EditStop
		a.JuncLen = 0
	}

	// If fragment matches the fully edited template entirely. Set junc len 0
	if T[0].All() {
		a.JuncEnd = a.JuncStart
		a.EditStop = a.JuncStart
		a.JuncLen = 0
		a.JuncSeq = ""
	}

	a.ReadCount = frag.ReadCount
	a.Norm = frag.Norm
	a.EditStop += int(tmpl.EditOffset)
	a.JuncStart += int(tmpl.EditOffset)
	a.JuncEnd += int(tmpl.EditOffset)

	return a
}

func readInt64(b []byte) int {
	n := (uint32(b[0]) << 24) |
		(uint32(b[1]) << 16) |
		(uint32(b[2]) << 8) |
		uint32(b[3])
	return int(int32(n))
}

func writeInt64(b []byte, n int) {
	b[0] = byte(uint64(n) >> 24)
	b[1] = byte(uint64(n) >> 16)
	b[2] = byte(uint64(n) >> 8)
	b[3] = byte(uint64(n))
}

func (a *Alignment) UnmarshalBinary(buf []byte) error {
	a.EditStop = readInt64(buf[0:4])
	a.JuncStart = readInt64(buf[4:8])
	a.JuncEnd = readInt64(buf[8:12])
	a.JuncLen = readInt64(buf[12:16])
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

	writeInt64(buf[0:4], a.EditStop)
	writeInt64(buf[4:8], a.JuncStart)
	writeInt64(buf[8:12], a.JuncEnd)
	writeInt64(buf[12:16], a.JuncLen)
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
	labels = append(labels, "RD")

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

	_, err = w.Write([]byte(fmt.Sprintf("JSS: %d\nESS: %d\nJES: %d\nJunc Len: %d\n", a.JuncStart, a.EditStop, a.JuncEnd, a.JuncLen)))
	if err != nil {
		return err
	}

	if a.AltEditing > 0 {
		_, err = w.Write([]byte(fmt.Sprintf("Alt: %d\n", a.AltEditing)))
		if err != nil {
			return err
		}
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
