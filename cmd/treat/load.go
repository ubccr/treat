// Copyright 2015 TREAT Authors. All rights reserved.
// Use of this source code is governed by a BSD style
// license that can be found in the LICENSE file.

package main

import (
    "os"
    "strings"
    "math"
    "fmt"
    "regexp"
    "bufio"
    "strconv"
    "path/filepath"
    "encoding/binary"
    "github.com/Sirupsen/logrus"
    "github.com/aebruno/gofasta"
    "github.com/ubccr/treat"
    "github.com/boltdb/bolt"
)

type LoadOptions struct {
    Gene          string
    EditBase      string
    TemplatePath  string
    FragmentPath  []string
    Primer5       string
    Primer3       string
    Replicate     int
    Norm          float64
    SkipFrags     bool
    ExcludeSnps   bool
    Fastx         bool
    Force         bool
    Collapse      bool
}

var readsPattern = regexp.MustCompile(`\s*merge_count=(\d+)\s*`)
var fastxPattern = regexp.MustCompile(`^\d+\-(\d+)$`)

func MergeCount(rec *gofasta.SeqRecord, options *LoadOptions) (uint32) {
    if options.Fastx {
        matches := fastxPattern.FindStringSubmatch(rec.Id)
        if len(matches) == 2 {
            count, err := strconv.Atoi(matches[1])
            if err != nil {
                logrus.Fatalf("Invalid FASTX header line found: %s", rec.Id)
            }

            return uint32(count)
        }
        logrus.Fatalf("Invalid FASTX header line found: %s", rec.Id)
    }

    mergeCount := 1
    matches := readsPattern.FindStringSubmatch(rec.Id)
    if len(matches) == 2 {
        count, err := strconv.Atoi(matches[1])
        if err == nil {
            mergeCount = count
        }
    }

    return uint32(mergeCount)
}

func TotalReads(path string, options *LoadOptions) (uint32) {
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
    for _,path := range(options.FragmentPath) {
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
            seqNoPrimer,err := frag.SequenceNoPrimer(tmpl.Primer5, tmpl.Primer3)
            if err != nil {
                logrus.Fatal(err)
            }
            mergeFrags[seqNoPrimer] = append(mergeFrags[seqNoPrimer], frag)
            mergeCounts[seqNoPrimer] += frag.ReadCount
        }

        outPath := strings.Replace(path, filepath.Ext(path), "", 1)
        outFile, err := os.Create(outPath+"-primer-collapsed.fa")
        if err != nil {
            logrus.Fatal(err)
        }

        defer outFile.Close()
        writer := bufio.NewWriter(outFile)

        i := 1
        for seq, frags := range mergeFrags {
            max := frags[0].ReadCount
            frag := frags[0]

            for _,f := range frags {
                if f.ReadCount > max {
                    max = f.ReadCount
                    frag = f
                }
            }

            writer.WriteString(">"+strconv.Itoa(i)+"-"+strconv.Itoa(int(mergeCounts[seq])))
            writer.WriteString("\n")
            writer.WriteString(frag.String())
            writer.WriteString("\n")
            i++
        }
        writer.Flush()
    }
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
        return
    }

    // Compute Edit Stop Site
    tmpl.EditStop = uint32((tmpl.Len()-1)-tmpl.Primer3)
    for j := int(tmpl.EditStop); j >= int(tmpl.Primer5); j-- {
        if tmpl.EditSite[0][j] != tmpl.EditSite[1][j] {
            tmpl.EditStop = uint32((tmpl.Len()-1)-j)
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
        for _,path := range(options.FragmentPath) {
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

    for _,path := range(options.FragmentPath) {
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
        sample = strings.Replace(sample, " ", "_", -1)

        akey := &treat.AlignmentKey{options.Gene, sample}
        key, err := akey.MarshalBinary()
        if err != nil {
            logrus.Fatal(err)
        }

        err = db.Update(func(tx *bolt.Tx) error {
            b := tx.Bucket([]byte(BUCKET_ALIGNMENTS))
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

            if count % 100 == 0 {
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

        if count % 100 != 0 {
            if err := tx.Commit(); err != nil {
                logrus.Fatal(err)
            }
        }

        logrus.Printf("Loaded %d fragment sequences for sample %s", count, sample)
    }
}
