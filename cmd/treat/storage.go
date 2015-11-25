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
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/Sirupsen/logrus"
	"github.com/aebruno/gofasta"
	"github.com/boltdb/bolt"
	"github.com/ubccr/treat"
)

const BUCKET_ALIGNMENTS = "alignments"
const BUCKET_TEMPLATES = "templates"
const BUCKET_FRAGMENTS = "fragments"
const BUCKET_META = "meta"
const STORAGE_VERSION_KEY = "version"
const STORAGE_VERSION = 0.2

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
	return openBolt(dbpath, 0644, &bolt.Options{Timeout: 1 * time.Second, ReadOnly: true})
}

func NewStorageWrite(dbpath string) (*Storage, error) {
	return openBolt(dbpath, 0644, &bolt.Options{Timeout: 1 * time.Second})
}

func openBolt(dbpath string, mode os.FileMode, options *bolt.Options) (*Storage, error) {
	db, err := bolt.Open(dbpath, mode, options)
	if err != nil {
		return nil, fmt.Errorf("Failed to open database %s - %s", dbpath, err)
	}

	storage := &Storage{DB: db}

	if options.ReadOnly {
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
	}

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
		if b == nil {
			return fmt.Errorf("database error. templates bucket does not exist!")
		}

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
		if b == nil {
			return fmt.Errorf("database error. templates bucket does not exist!")
		}

		v := b.Get([]byte(gene))
		if v == nil {
			return fmt.Errorf("database error. template not found for gene: %s", gene)
		}

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

func (s *Storage) Initialize() error {
	err := s.DB.Update(func(tx *bolt.Tx) error {
		_, err := tx.CreateBucketIfNotExists([]byte(BUCKET_ALIGNMENTS))
		if err != nil {
			return err
		}

		_, err = tx.CreateBucketIfNotExists([]byte(BUCKET_FRAGMENTS))
		if err != nil {
			return err
		}

		b, err := tx.CreateBucketIfNotExists([]byte(BUCKET_META))
		if err != nil {
			return err
		}

		versionBytes := b.Get([]byte(STORAGE_VERSION_KEY))
		if versionBytes != nil {
			version := math.Float64frombits(binary.BigEndian.Uint64(versionBytes))
			if version != STORAGE_VERSION {
				return fmt.Errorf("treat version mismatch. Must re-load data using same version of treat. %.1f != %.1f", version, STORAGE_VERSION)
			}
		}

		vbuf := make([]byte, 8)
		binary.BigEndian.PutUint64(vbuf, math.Float64bits(STORAGE_VERSION))
		err = b.Put([]byte(STORAGE_VERSION_KEY), vbuf)
		if err != nil {
			return err
		}

		b, err = tx.CreateBucketIfNotExists([]byte(BUCKET_TEMPLATES))
		if err != nil {
			return err
		}

		return nil
	})

	return err
}

func (s *Storage) InitSample(akey *treat.AlignmentKey, force bool) error {
	key, err := akey.MarshalBinary()
	if err != nil {
		return err
	}

	err = s.DB.Update(func(tx *bolt.Tx) error {
		b := tx.Bucket([]byte(BUCKET_ALIGNMENTS))
		if b == nil {
			return fmt.Errorf("database error. alignments bucket does not exist!")
		}

		_, err := b.CreateBucket(key)
		if err != nil {
			if !force {
				return fmt.Errorf("Data already exists for gene %s and sample %s. Use --force to force delete data and reload (error: %s)", akey.Gene, akey.Sample, err)
			}

			logrus.WithFields(logrus.Fields{
				"gene":   akey.Gene,
				"sample": akey.Sample,
			}).Warn("Deleting existing alignment data")
			err = b.DeleteBucket(key)
			if err != nil {
				return fmt.Errorf("database error. failed to delete nested alignment bucket: %s", err)
			}

			_, err := b.CreateBucket(key)
			if err != nil {
				return fmt.Errorf("database error. failed to create nested alignment bucket: %s", err)
			}
		}

		b = tx.Bucket([]byte(BUCKET_FRAGMENTS))
		if b == nil {
			return fmt.Errorf("database error. fragments bucket does not exist!")
		}

		_, err = b.CreateBucket(key)
		if err != nil {
			if !force {
				return fmt.Errorf("Data already exists for gene %s and sample %s. Please delete database and reload (error: %s)", akey.Gene, akey.Sample, err)
			}

			logrus.WithFields(logrus.Fields{
				"gene":   akey.Gene,
				"sample": akey.Sample,
			}).Warn("Deleting existing fragment data")
			err = b.DeleteBucket(key)
			if err != nil {
				return fmt.Errorf("database error. failed to delete nested fragment bucket: %s", err)
			}

			_, err := b.CreateBucket(key)
			if err != nil {
				return fmt.Errorf("database error. failed to create nested fragment bucket: %s", err)
			}
		}

		return nil
	})

	return err
}

func (s *Storage) ImportSample(path string, options *LoadOptions) (*treat.AlignmentKey, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	tmpl, err := s.GetTemplate(options.Gene)
	if err != nil {
		return nil, err
	}

	fname := filepath.Base(path)
	sample := fname[:len(fname)-len(filepath.Ext(path))]
	// clean up sample name
	sample = strings.Replace(sample, "-primer-collapsed", "", 1)
	sample = strings.Replace(sample, " ", "_", -1)

	akey := &treat.AlignmentKey{options.Gene, sample}
	key, err := akey.MarshalBinary()
	if err != nil {
		return nil, err
	}

	err = s.InitSample(akey, options.Force)
	if err != nil {
		return nil, err
	}

	var tx *bolt.Tx
	var alnBucket *bolt.Bucket
	var fragBucket *bolt.Bucket
	count := 0

	logrus.Printf("Processing fragments for sample name: %s", sample)
	if options.SkipFrags {
		logrus.Info("not storing raw fragment reads")
	}

	for rec := range gofasta.SimpleParser(f) {

		if count%100 == 0 {
			if count > 0 {
				if err := tx.Commit(); err != nil {
					return nil, err
				}
			}
			tx, err = s.DB.Begin(true)
			if err != nil {
				return nil, err
			}
			alnBucket = tx.Bucket([]byte(BUCKET_ALIGNMENTS)).Bucket(key)
			fragBucket = tx.Bucket([]byte(BUCKET_FRAGMENTS)).Bucket(key)
		}

		frag := treat.NewFragment(rec.Id, rec.Seq, treat.FORWARD, rune(options.EditBase[0]))
		aln := treat.NewAlignment(frag, tmpl, options.ExcludeSnps)

		id, _ := alnBucket.NextSequence()
		kbytes := make([]byte, 8)
		binary.BigEndian.PutUint64(kbytes, id)

		data, err := aln.MarshalBinary()
		if err != nil {
			return nil, err
		}

		err = alnBucket.Put(kbytes, data)
		if err != nil {
			return nil, err
		}

		if !options.SkipFrags {
			data, err = frag.MarshalBytes()
			if err != nil {
				return nil, err
			}

			err = fragBucket.Put(kbytes, data)
			if err != nil {
				return nil, err
			}
		}

		count++
	}

	// final transaction commit
	if err := tx.Commit(); err != nil {
		return nil, err
	}

	logrus.Printf("Loaded %d fragment sequences for sample %s", count, sample)

	return akey, nil
}

func (s *Storage) NormalizeSample(akey *treat.AlignmentKey, norm float64) error {
	key, err := akey.MarshalBinary()
	if err != nil {
		return err
	}

	total := 0
	err = s.DB.View(func(tx *bolt.Tx) error {
		ab := tx.Bucket([]byte(BUCKET_ALIGNMENTS))
		if ab == nil {
			return fmt.Errorf("database error. alignments bucket does not exist!")
		}

		b := ab.Bucket(key)
		if ab == nil {
			return fmt.Errorf("database error. key not found in alignments bucket")
		}

		c := b.Cursor()
		for ak, av := c.First(); ak != nil; ak, av = c.Next() {
			a := new(treat.Alignment)
			a.Id = binary.BigEndian.Uint64(ak)
			err := a.UnmarshalBinary(av)
			if err != nil {
				return err
			}

			// Only count Standard Reads
			if a.HasMutation == 0 {
				total += int(a.ReadCount)
			}
		}

		return nil
	})

	if err != nil {
		return err
	}

	scale := 1.0
	if total > 0 {
		scale = norm / float64(total)
	}

	logrus.Printf("Processing sample %s using normalized scaling factor: %.4f", akey.Sample, scale)
	err = s.DB.Update(func(tx *bolt.Tx) error {
		ab := tx.Bucket([]byte(BUCKET_ALIGNMENTS))
		b := ab.Bucket(key)

		c := b.Cursor()
		for ak, av := c.First(); ak != nil; ak, av = c.Next() {
			a := new(treat.Alignment)
			a.Id = binary.BigEndian.Uint64(ak)
			err := a.UnmarshalBinary(av)
			if err != nil {
				return err
			}

			if a.HasMutation == uint8(0) {
				a.Norm = scale * float64(a.ReadCount)
			}

			data, err := a.MarshalBinary()
			if err != nil {
				return err
			}

			err = b.Put(ak, data)
			if err != nil {
				return err
			}
		}

		return nil
	})

	if err != nil {
		return err
	}

	return nil
}
