package main

import (
    "fmt"
    "os"
    "log"
    "regexp"
    "strconv"
    "path/filepath"
    "github.com/aebruno/gofasta"
    "github.com/ubccr/treat"
    "github.com/boltdb/bolt"
)

type LoadOptions struct {
    Gene          string
    TemplatePath  string
    FragmentPath  []string
    Primer5       int
    Primer3       int
    Replicate     int
    Norm          float64
    GrnaPath      string
}

var readsPattern = regexp.MustCompile(`\s*merge_count=(\d+)\s*`)

func Load(dbpath string, options *LoadOptions) {
    if len(options.Gene) == 0 {
        log.Fatalln("ERROR Gene name is required")
    }
    if len(options.TemplatePath) == 0 {
        log.Fatalln("ERROR Please provide path to templates file")
    }
    if options.FragmentPath == nil || len(options.FragmentPath) == 0 {
        log.Fatalln("ERROR Please provide 1 or more fragment files to load")
    }

    tmpl, err := treat.NewTemplateFromFasta(options.TemplatePath, treat.FORWARD, 't')
    if err != nil {
        log.Fatalln(err)
    }

    db, err := bolt.Open(dbpath, 0644, nil)
    if err != nil {
        log.Fatal(err)
    }

    err = db.Update(func(tx *bolt.Tx) error {
        _, err := tx.CreateBucketIfNotExists([]byte("alignments"))
        if err != nil {
            return err
        }
        return nil
    })

    if err != nil {
        log.Fatal(err)
    }

    for _,path := range(options.FragmentPath) {
        f, err := os.Open(path)
        if err != nil {
            log.Fatal(err)
        }

        fname := filepath.Base(path)
        sample := fname[:len(fname)-len(filepath.Ext(path))]

        var tx *bolt.Tx
        count := 0

        for rec := range gofasta.SimpleParser(f) {

            if count % 100 == 0 {
                if count > 0 {
                    if err := tx.Commit(); err != nil {
                        log.Fatal(err)
                    }
                }
                tx, err = db.Begin(true)
                if err != nil {
                    log.Fatal(err)
                }
            }

            mergeCount := 0
            matches := readsPattern.FindStringSubmatch(rec.Id)
            if len(matches) == 2 {
                mergeCount, err = strconv.Atoi(matches[1])
                if err != nil {
                    mergeCount = 0
                }
            }

            frag := treat.NewFragment(rec.Id, rec.Seq, treat.FORWARD, treat.ReadCountType(mergeCount), 't')
            aln := treat.NewAlignment(frag, tmpl, options.Primer5, options.Primer3)


            bucket := tx.Bucket([]byte("alignments"))
            id, _ := bucket.NextSequence()

            key := fmt.Sprintf("%s;%s;%d;%d", options.Gene, sample, options.Replicate, id)

            data, err := aln.Bytes()
            if err != nil {
                log.Fatal(err)
            }

            err = bucket.Put([]byte(key), data)
            if err != nil {
                log.Fatal(err)
            }

            count++
        }

        if count % 100 != 0 {
            if err := tx.Commit(); err != nil {
                log.Fatal(err)
            }
        }
    }
}
