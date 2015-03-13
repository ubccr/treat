package main

import (
    "log"
    "fmt"
    "bytes"
    "strings"
    "strconv"
    "math"
    "github.com/boltdb/bolt"
    "github.com/ubccr/treat"
)

type SearchFields struct {
    Gene          string
    Sample        string
    Replicate     int
    EditStop      int
    JuncEnd       int
    Offset        int
    Limit         int
    HasMutation   bool
    HasAlt        bool
    All           bool
    GrnaEdit      []int
    GrnaJunc      []int
}

func (fields *SearchFields) HasKeyMatch(k string) bool {
    key := strings.Split(string(k), ";")
    replicate,_ := strconv.Atoi(key[2])

    if fields.Replicate > 0 && replicate != fields.Replicate {
        return false
    }

    if len(fields.Sample) > 0 && fields.Sample != key[1] {
        return false
    }
    if fields.Replicate > 0 && fields.Replicate != replicate {
        return false
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

    if fields.EditStop > 0 && uint64(fields.EditStop) != a.EditStop {
        return false
    }
    if fields.JuncEnd > 0 && uint64(fields.JuncEnd) != a.JuncEnd {
        return false
    }
    if fields.HasAlt && a.AltEditing == -1 {
        return false
    }
    gflag := false
    for _, g := range(fields.GrnaEdit) {
        if a.GrnaEdit.Bit(g) == 0 {
            gflag = true
        }
    }
    for _, g := range(fields.GrnaJunc) {
        if a.GrnaJunc.Bit(g) == 0 {
            gflag = true
        }
    }
    if gflag {
        return false
    }

    return true
}

// From: https://gist.github.com/DavidVaini/10308388
func Round(f float64) float64 {
    return math.Floor(f + .5)
}

func RoundPlus(f float64, places int) (float64) {
    shift := math.Pow(10, float64(places))
    return Round(f * shift) / shift;
}

func PrintResult(k string, a *treat.Alignment) {
    key := strings.Split(string(k), ";")

    fmt.Println(strings.Join([]string{
        key[0],
        key[1],
        key[2],
        key[3],
        fmt.Sprintf("%.4f", RoundPlus(a.Norm, 4)),
        fmt.Sprintf("%d", a.ReadCount),
        fmt.Sprintf("%d", a.AltEditing),
        fmt.Sprintf("%d", a.HasMutation),
        fmt.Sprintf("%d", a.EditStop),
        fmt.Sprintf("%d", a.JuncEnd),
        fmt.Sprintf("%d", a.JuncLen),
        a.GrnaEditString()}, "\t"))
}

func Search(dbpath string, fields *SearchFields) {
    db, err := bolt.Open(dbpath, 0644, nil)
    if err != nil {
        log.Fatal(err)
    }

    count := 0
    offset := 0

    db.View(func(tx *bolt.Tx) error {
        c := tx.Bucket([]byte("alignments")).Cursor()

        prefix := ""
        if len(fields.Gene) > 0 {
            prefix +=  fields.Gene
            if len(fields.Sample) > 0 {
                prefix += ";"+fields.Sample
                if fields.Replicate > 0 {
                    prefix =  fmt.Sprintf("%s;%d", prefix, fields.Replicate)
                }
            }
        }

        if len(prefix) > 0 {
            for k, v := c.Seek([]byte(prefix)); bytes.HasPrefix(k, []byte(prefix)); k, v = c.Next() {
                key := string(k)
                if !fields.HasKeyMatch(key) {
                    continue
                }
                a, err := treat.NewAlignmentFromBytes(v)
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

                PrintResult(key, a)
                count++
                offset++
            }
        } else {
            db.View(func(tx *bolt.Tx) error {
                b := tx.Bucket([]byte("alignments"))
                b.ForEach(func(k, v []byte) error {
                    key := string(k)
                    if !fields.HasKeyMatch(key) {
                        return nil
                    }

                    a, err := treat.NewAlignmentFromBytes(v)
                    if err != nil {
                        return err
                    }

                    if !fields.HasMatch(a) {
                        return nil
                    }

                    if fields.Offset > 0 && offset < fields.Offset {
                        offset++
                        return nil
                    }

                    if fields.Limit > 0 && count >= fields.Limit {
                        return nil
                    }

                    PrintResult(key, a)
                    count++
                    offset++

                    return nil
                })
                return nil
            })
        }

        return nil
    })
}
