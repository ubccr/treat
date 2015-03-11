package main

import (
    "log"
    "fmt"
    "bytes"
    "strings"
    "strconv"
    "github.com/boltdb/bolt"
    "github.com/ubccr/treat"
)

type SearchFields struct {
    Gene          string
    Sample        string
    Replicate     int
    EditStop      int
    JuncEnd       int
    HasMutation   bool
}

func Search(dbpath string, fields *SearchFields) {
    db, err := bolt.Open(dbpath, 0644, nil)
    if err != nil {
        log.Fatal(err)
    }

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
                key := strings.Split(string(k), ";")
                replicate,_ := strconv.Atoi(key[2])

                if fields.Replicate > 0 && replicate != fields.Replicate {
                    continue
                }

                a, err := treat.NewAlignmentFromBytes(v)
                if err != nil {
                    return err
                }

                if fields.HasMutation && a.HasMutation == 0 {
                    continue
                } else if !fields.HasMutation && a.HasMutation == 1 {
                    continue
                }

                if fields.EditStop > 0 && uint64(fields.EditStop) != a.EditStop {
                    continue
                }
                if fields.JuncEnd > 0 && uint64(fields.JuncEnd) != a.JuncEnd {
                    continue
                }

                fmt.Println(strings.Join([]string{key[0],key[1],key[2],key[3],fmt.Sprintf("%d", a.EditStop),fmt.Sprintf("%d", a.JuncEnd)}, "\t"))
            }
        } else {
            db.View(func(tx *bolt.Tx) error {
                b := tx.Bucket([]byte("alignments"))
                b.ForEach(func(k, v []byte) error {
                    key := strings.Split(string(k), ";")
                    replicate,_ := strconv.Atoi(key[2])
                    if len(fields.Sample) > 0 && fields.Sample != key[1] {
                        return nil
                    }
                    if fields.Replicate > 0 && fields.Replicate != replicate {
                        return nil
                    }


                    a, err := treat.NewAlignmentFromBytes(v)
                    if err != nil {
                        return err
                    }

                    if fields.EditStop > 0 && uint64(fields.EditStop) != a.EditStop {
                        return nil
                    }
                    if fields.JuncEnd > 0 && uint64(fields.JuncEnd) != a.JuncEnd {
                        return nil
                    }

                    if fields.HasMutation && a.HasMutation == 0 {
                        return nil
                    } else if a.HasMutation == 1 {
                        return nil
                    }

                    fmt.Println(strings.Join([]string{key[0],key[1],key[2],key[3],fmt.Sprintf("%d", a.EditStop),fmt.Sprintf("%d", a.JuncEnd)}, "\t"))

                    return nil
                })
                return nil
            })
        }

        return nil
    })
}
