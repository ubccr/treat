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

package main

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"math"
	"time"

	"github.com/boltdb/bolt"
	"github.com/ubccr/treat"
)

const BUCKET_ALIGNMENTS = "alignments"
const BUCKET_TEMPLATES = "templates"
const BUCKET_FRAGMENTS = "fragments"
const BUCKET_META = "meta"
const STORAGE_VERSION_KEY = "version"
const STORAGE_VERSION = 0.1

type Storage struct {
	DB      *bolt.DB
	version float64
}

type SearchFields struct {
	Gene        string   `schema:"gene"`
	Sample      []string `schema:"sample"`
	EditStop    int      `schema:"edit_stop"`
	JuncEnd     int      `schema:"junc_end"`
	JuncLen     int      `schema:"junc_len"`
	Offset      int      `schema:"offset"`
	Limit       int      `schema:"limit"`
	HasMutation bool     `schema:"has_mutation"`
	HasAlt      bool     `schema:"has_alt"`
	All         bool     `schema:"all"`
	AltRegion   int      `schema:"alt"`
}

type AlignmentResults []*treat.Alignment

func (s AlignmentResults) Len() int      { return len(s) }
func (s AlignmentResults) Swap(i, j int) { s[i], s[j] = s[j], s[i] }

type ByReadCount struct{ AlignmentResults }

func (s ByReadCount) Less(i, j int) bool {
	return s.AlignmentResults[i].ReadCount > s.AlignmentResults[j].ReadCount
}

func (fields *SearchFields) HasSample(val string) bool {
	for _, s := range fields.Sample {
		if s == val {
			return true
		}
	}

	return false
}

func (fields *SearchFields) HasKeyMatch(k *treat.AlignmentKey) bool {
	if len(fields.Gene) > 0 && fields.Gene != k.Gene {
		return false
	}

	if len(fields.Sample) > 0 {
		match := false
		for _, s := range fields.Sample {
			if s == k.Sample {
				match = true
				break
			}
		}
		if !match {
			return false
		}
	}

	return true
}

func (fields *SearchFields) HasMatch(a *treat.Alignment) bool {
	if !fields.All {
		if fields.HasMutation && a.HasMutation == 0 {
			return false
		} else if !fields.HasMutation && a.HasMutation == 1 {
			return false
		}
	}

	if fields.EditStop >= 0 && uint32(fields.EditStop) != a.EditStop {
		return false
	}
	if fields.JuncLen >= 0 && uint32(fields.JuncLen) != a.JuncLen {
		return false
	}
	if fields.JuncEnd >= 0 && uint32(fields.JuncEnd) != a.JuncEnd {
		return false
	}
	if fields.HasAlt && a.AltEditing == 0 {
		return false
	}
	if fields.AltRegion > 0 && uint8(fields.AltRegion) != a.AltEditing {
		return false
	}

	return true
}

// From: https://gist.github.com/DavidVaini/10308388
func Round(f float64) float64 {
	return math.Floor(f + .5)
}

func RoundPlus(f float64, places int) float64 {
	shift := math.Pow(10, float64(places))
	return Round(f*shift) / shift
}

func NewStorage(dbpath string) (*Storage, error) {
	db, err := bolt.Open(dbpath, 0600, &bolt.Options{Timeout: 1 * time.Second})
	if err != nil {
		return nil, fmt.Errorf("Failed to open database %s - %s", dbpath, err)
	}

	storage := &Storage{DB: db}

	err = db.View(func(tx *bolt.Tx) error {
		b := tx.Bucket([]byte(BUCKET_META))
		if b == nil {
			return fmt.Errorf("Invalid db file. missing treat metadata")
		}

		versionBytes := b.Get([]byte(STORAGE_VERSION_KEY))
		if versionBytes == nil {
			return fmt.Errorf("Invalid db file. missing treat version")
		}
		storage.version = math.Float64frombits(binary.BigEndian.Uint64(versionBytes))

		return nil
	})

	if err != nil {
		return nil, err
	}

	return storage, nil
}

func (s *Storage) Search(fields *SearchFields, f func(k *treat.AlignmentKey, a *treat.Alignment)) error {
	count := 0
	offset := 0

	err := s.DB.View(func(tx *bolt.Tx) error {
		b := tx.Bucket([]byte(BUCKET_ALIGNMENTS))
		c := b.Cursor()

		for k, _ := c.First(); k != nil; k, _ = c.Next() {
			key := new(treat.AlignmentKey)
			key.UnmarshalBinary(k)
			if !fields.HasKeyMatch(key) {
				continue
			}

			bucket := c.Bucket().Bucket(k).Cursor()

			for ak, av := bucket.First(); ak != nil; ak, av = bucket.Next() {
				a := new(treat.Alignment)
				a.Id = binary.BigEndian.Uint64(ak)
				err := a.UnmarshalBinary(av)
				if err != nil {
					return err
				}

				if !fields.HasMatch(a) {
					continue
				}

				if fields.Offset > 0 && offset < fields.Offset {
					offset++
					continue
				}

				if fields.Limit > 0 && count >= fields.Limit {
					return nil
				}

				f(key, a)
				count++
				offset++
			}
		}

		return nil
	})

	return err
}

func (s *Storage) PutTemplate(gene string, tmpl *treat.Template) error {
	err := s.DB.Update(func(tx *bolt.Tx) error {
		b := tx.Bucket([]byte(BUCKET_TEMPLATES))

		data, err := tmpl.MarshalBytes()
		if err != nil {
			return err
		}

		err = b.Put([]byte(gene), data)
		return err
	})

	return err
}

func (s *Storage) GetTemplate(gene string) (*treat.Template, error) {
	var tmpl *treat.Template
	err := s.DB.View(func(tx *bolt.Tx) error {
		b := tx.Bucket([]byte(BUCKET_TEMPLATES))
		v := b.Get([]byte(gene))

		t := new(treat.Template)
		err := t.UnmarshalBytes(v)
		if err != nil {
			return err
		}

		tmpl = t

		return nil
	})

	if err != nil {
		return nil, err
	}

	return tmpl, nil
}

func (s *Storage) TemplateMap() (map[string]*treat.Template, error) {
	templates := make(map[string]*treat.Template, 0)
	err := s.DB.View(func(tx *bolt.Tx) error {
		b := tx.Bucket([]byte(BUCKET_TEMPLATES))

		b.ForEach(func(k, v []byte) error {

			tmpl := new(treat.Template)
			err := tmpl.UnmarshalBytes(v)
			if err != nil {
				return err
			}

			templates[string(k)] = tmpl

			return nil
		})

		return nil
	})

	if err != nil {
		return nil, err
	}

	return templates, nil
}

func (s *Storage) Genes() ([]string, error) {
	genes := make([]string, 0)
	err := s.DB.View(func(tx *bolt.Tx) error {
		b := tx.Bucket([]byte(BUCKET_TEMPLATES))

		b.ForEach(func(k, v []byte) error {
			genes = append(genes, string(k))
			return nil
		})

		return nil
	})

	if err != nil {
		return nil, err
	}

	return genes, nil
}

func (s *Storage) Samples(gene string) ([]string, error) {
	samples := make([]string, 0)
	gbytes := []byte(gene)

	err := s.DB.View(func(tx *bolt.Tx) error {
		c := tx.Bucket([]byte(BUCKET_ALIGNMENTS)).Cursor()
		for k, _ := c.Seek(gbytes); bytes.HasPrefix(k, gbytes); k, _ = c.Next() {
			key := new(treat.AlignmentKey)
			key.UnmarshalBinary(k)
			samples = append(samples, key.Sample)
		}

		return nil
	})

	if err != nil {
		return nil, err
	}

	return samples, nil
}

func (s *Storage) GetAlignment(k *treat.AlignmentKey, id uint64) (*treat.Alignment, error) {
	key, err := k.MarshalBinary()
	if err != nil {
		return nil, err
	}

	buf := make([]byte, 8)
	binary.BigEndian.PutUint64(buf, id)

	var alignment *treat.Alignment
	err = s.DB.View(func(tx *bolt.Tx) error {
		b := tx.Bucket([]byte(BUCKET_ALIGNMENTS)).Bucket(key)
		v := b.Get(buf)
		if v != nil {
			a := new(treat.Alignment)
			err := a.UnmarshalBinary(v)
			if err != nil {
				return err
			}

			alignment = a
		}

		return nil
	})

	if err != nil {
		return nil, err
	}

	return alignment, nil
}

func (s *Storage) GetFragment(k *treat.AlignmentKey, id uint64) (*treat.Fragment, error) {
	key, err := k.MarshalBinary()
	if err != nil {
		return nil, err
	}

	buf := make([]byte, 8)
	binary.BigEndian.PutUint64(buf, id)

	var frag *treat.Fragment
	err = s.DB.View(func(tx *bolt.Tx) error {
		b := tx.Bucket([]byte(BUCKET_FRAGMENTS)).Bucket(key)
		v := b.Get(buf)
		if v != nil {
			f := new(treat.Fragment)
			err := f.UnmarshalBytes(v)
			if err != nil {
				return err
			}

			frag = f
		}

		return nil
	})

	if err != nil {
		return nil, err
	}

	return frag, nil
}
