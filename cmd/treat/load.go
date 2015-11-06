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
	"bufio"
	"encoding/binary"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"

	"github.com/Sirupsen/logrus"
	"github.com/aebruno/gofasta"
	"github.com/boltdb/bolt"
	"github.com/ubccr/treat"
)

type LoadOptions struct {
	Gene         string
	EditBase     string
	TemplatePath string
	FragmentPath []string
	Primer5      string
	Primer3      string
	Replicate    int
	Norm         float64
	SkipFrags    bool
	ExcludeSnps  bool
	Force        bool
	Collapse     bool
}

// fastx_collapser will output fasta headers in the format [id]-[count]
//  for example: 
//      > 132-2082
//
// This regex will parse any fasta header (including those from
// fastx_collapser) where the collapsed read count is found at the end of the
// string separated by either a '-' or an '_'. For example, all these 
// are acceptable and the collapsed count would = 2082:
//     > SAMPLE1_GENE_123432_2082
//     > 132-2082
//     > GENE-88772-2082
var fastxPattern = regexp.MustCompile(`.+[_\-](\d+)$`)

func MergeCount(rec *gofasta.SeqRecord, options *LoadOptions) uint32 {
	// First try and parse fastx-collapser header lines
	matches := fastxPattern.FindStringSubmatch(rec.Id)
	if len(matches) == 2 {
		count, err := strconv.Atoi(matches[1])
		if err != nil {
			logrus.WithFields(logrus.Fields{
				"error": err,
				"id":    rec.Id,
			}).Warn("Looks like a FASTX header line but failed to parse. Skipping.")
			return uint32(1)
		}

		return uint32(count)
	}

	return uint32(1)
}

func TotalReads(path string, options *LoadOptions) uint32 {
	total := uint32(0)

	f, err := os.Open(path)
	if err != nil {
		logrus.Fatal(err)
	}
	defer f.Close()

	for rec := range gofasta.SimpleParser(f) {
		total += MergeCount(rec, options)
	}

	return total
}

func collapse(tmpl *treat.Template, options *LoadOptions) {
	logrus.Info("Collapsing fragment files excluding primer regions")
	collapsedFiles := make([]string, 0)
	for _, path := range options.FragmentPath {
		logrus.Printf("Processing file: %s", path)
		inFile, err := os.Open(path)
		if err != nil {
			logrus.Fatal(err)
		}
		defer inFile.Close()

		mergeFrags := make(map[string][]*treat.Fragment)
		mergeCounts := make(map[string]uint32)

		for rec := range gofasta.SimpleParser(inFile) {
			mergeCount := MergeCount(rec, options)
			frag := treat.NewFragment(rec.Id, rec.Seq, treat.FORWARD, mergeCount, 0, rune(options.EditBase[0]))
			seqNoPrimer, err := frag.StripPrimer(tmpl.Primer5, tmpl.Primer3)
			if err != nil {
				logrus.Fatal(err)
			}
			mergeFrags[seqNoPrimer] = append(mergeFrags[seqNoPrimer], frag)
			mergeCounts[seqNoPrimer] += frag.ReadCount
		}

		outPath := path[:len(path)-len(filepath.Ext(path))]
		outPath += "-primer-collapsed.fasta"
		outFile, err := os.Create(outPath)
		if err != nil {
			logrus.Fatal(err)
		}

		defer outFile.Close()
		writer := bufio.NewWriter(outFile)

		i := 1
		for seq, frags := range mergeFrags {
			max := frags[0].ReadCount
			frag := frags[0]

			for _, f := range frags {
				if f.ReadCount > max {
					max = f.ReadCount
					frag = f
				}
			}

			writer.WriteString(">" + strconv.Itoa(i) + "-" + strconv.Itoa(int(mergeCounts[seq])))
			writer.WriteString("\n")
			writer.WriteString(frag.String())
			writer.WriteString("\n")
			i++
		}
		writer.Flush()
		collapsedFiles = append(collapsedFiles, outPath)
	}

	// Process new collapsed files
	options.FragmentPath = collapsedFiles
}

func Load(dbpath string, options *LoadOptions) {
	if len(options.Gene) == 0 {
		logrus.Fatal("Gene name is required")
	}
	if len(options.TemplatePath) == 0 {
		logrus.Fatal("Please provide path to templates file")
	}
	if options.FragmentPath == nil || len(options.FragmentPath) == 0 {
		logrus.Fatal("Please provide 1 or more fragment files to load")
	}
	if len(options.EditBase) != 1 {
		logrus.Fatal("Please provide the edit base")
	}

	options.Gene = strings.Replace(options.Gene, " ", "_", -1)

	tmpl, err := treat.NewTemplateFromFasta(options.TemplatePath, treat.FORWARD, rune(options.EditBase[0]))
	if err != nil {
		logrus.Fatalln(err)
	}

	if len(options.Primer3) > 0 {
		err = tmpl.SetPrimer3(options.Primer3)
		if err != nil {
			logrus.Fatalln(err)
		}
	}

	if len(options.Primer5) > 0 {
		err = tmpl.SetPrimer5(options.Primer5)
		if err != nil {
			logrus.Fatalln(err)
		}
	}

	if options.Collapse {
		collapse(tmpl, options)
	}

	// Compute Edit Stop Site
	tmpl.EditStop = uint32((tmpl.Len() - 1) - tmpl.Primer3)
	for j := int(tmpl.EditStop); j >= int(tmpl.Primer5); j-- {
		if tmpl.EditSite[0][j] != tmpl.EditSite[1][j] {
			tmpl.EditStop = uint32((tmpl.Len() - 1) - j)
			break
		}
	}
	if tmpl.EditStop > 0 {
		tmpl.EditStop--
	}
	logrus.Printf("Using template Edit Stop Site: %d", tmpl.EditStop)

	// Normalize read counts
	if options.Norm == 0 {
		total := uint32(0)
		for _, path := range options.FragmentPath {
			total += TotalReads(path, options)
		}
		options.Norm = float64(total) / float64(len(options.FragmentPath))
		logrus.Printf("Total reads across all samples: %d", total)
		logrus.Printf("Normalizing to average read count:: %.4f", options.Norm)
	}

	db, err := bolt.Open(dbpath, 0644, nil)
	if err != nil {
		logrus.Fatal(err)
	}

	err = db.Update(func(tx *bolt.Tx) error {
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

		data, err := tmpl.MarshalBytes()
		if err != nil {
			return err
		}

		err = b.Put([]byte(options.Gene), data)
		return err
	})

	if err != nil {
		logrus.Fatal(err)
	}

	for _, path := range options.FragmentPath {
		logrus.Printf("Computing total read count for file: %s", path)
		total := TotalReads(path, options)
		scale := options.Norm / float64(total)
		logrus.Printf("Total reads for file: %d", total)
		logrus.Printf("Normalized scaling factor: %.4f", scale)

		f, err := os.Open(path)
		if err != nil {
			logrus.Fatal(err)
		}
		defer f.Close()

		fname := filepath.Base(path)
		sample := fname[:len(fname)-len(filepath.Ext(path))]
		// clean up sample name
		sample = strings.Replace(sample, "-primer-collapsed", "", 1)
		sample = strings.Replace(sample, " ", "_", -1)

		akey := &treat.AlignmentKey{options.Gene, sample}
		key, err := akey.MarshalBinary()
		if err != nil {
			logrus.Fatal(err)
		}

		err = db.Update(func(tx *bolt.Tx) error {
			b := tx.Bucket([]byte(BUCKET_ALIGNMENTS))
			if b == nil {
				return fmt.Errorf("database error. alignments bucket does not exist!")
			}

			_, err := b.CreateBucket(key)
			if err != nil {
				if !options.Force {
					return fmt.Errorf("Data already exists for gene %s and sample %s. Use --force to force delete data and reload (error: %s)", akey.Gene, akey.Sample, err)
				}

				logrus.Printf("Deleting existing alignment data for gene %s and sample %s", akey.Gene, akey.Sample)
				err = b.DeleteBucket(key)
				if err != nil {
					return fmt.Errorf("database error. failed to delete bucket: %s", err)
				}

				_, err := b.CreateBucket(key)
				if err != nil {
					return fmt.Errorf("database error. failed to create bucket: %s", err)
				}
			}

			b = tx.Bucket([]byte(BUCKET_FRAGMENTS))
			if b == nil {
				return fmt.Errorf("database error. fragments bucket does not exist!")
			}

			_, err = b.CreateBucket(key)
			if err != nil {
				if !options.Force {
					return fmt.Errorf("Data already exists for gene %s and sample %s. Please delete database and reload (error: %s)", akey.Gene, akey.Sample, err)
				}

				logrus.Printf("Deleting existing fragment data for gene %s and sample %s", akey.Gene, akey.Sample)
				err = b.DeleteBucket(key)
				if err != nil {
					return fmt.Errorf("database error. failed to delete bucket: %s", err)
				}

				_, err := b.CreateBucket(key)
				if err != nil {
					return fmt.Errorf("database error. failed to create bucket: %s", err)
				}
			}

			return nil
		})

		if err != nil {
			logrus.Fatalf("%s", err)
		}

		var tx *bolt.Tx
		var alnBucket *bolt.Bucket
		var fragBucket *bolt.Bucket
		count := 0

		logrus.Printf("Processing fragments for sample name : %s", sample)
		if options.SkipFrags {
			logrus.Print("[info] not storing raw fragment reads")
		}

		for rec := range gofasta.SimpleParser(f) {

			if count%100 == 0 {
				if count > 0 {
					if err := tx.Commit(); err != nil {
						logrus.Fatal(err)
					}
				}
				tx, err = db.Begin(true)
				if err != nil {
					logrus.Fatal(err)
				}
				alnBucket = tx.Bucket([]byte(BUCKET_ALIGNMENTS)).Bucket(key)
				fragBucket = tx.Bucket([]byte(BUCKET_FRAGMENTS)).Bucket(key)
			}

			mergeCount := MergeCount(rec, options)
			norm := scale * float64(mergeCount)

			frag := treat.NewFragment(rec.Id, rec.Seq, treat.FORWARD, mergeCount, norm, rune(options.EditBase[0]))
			aln := treat.NewAlignment(frag, tmpl, options.ExcludeSnps)

			id, _ := alnBucket.NextSequence()
			kbytes := make([]byte, 8)
			binary.BigEndian.PutUint64(kbytes, id)

			data, err := aln.MarshalBinary()
			if err != nil {
				logrus.Fatal(err)
			}

			err = alnBucket.Put(kbytes, data)
			if err != nil {
				logrus.Fatal(err)
			}

			if !options.SkipFrags {
				data, err = frag.MarshalBytes()
				if err != nil {
					logrus.Fatal(err)
				}

				err = fragBucket.Put(kbytes, data)
				if err != nil {
					logrus.Fatal(err)
				}
			}

			count++
		}

		// final transaction commit
		if err := tx.Commit(); err != nil {
			logrus.Fatal(err)
		}

		logrus.Printf("Loaded %d fragment sequences for sample %s", count, sample)
	}
}
