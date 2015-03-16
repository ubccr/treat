package treat

import (
    "fmt"
    "time"
    "bytes"
    "strings"
    "strconv"
    "math"
    "github.com/boltdb/bolt"
)

type Storage struct {
    DB         *bolt.DB
}

type SearchFields struct {
    Gene          string      `schema:"gene"`
    Sample        []string    `schema:"sample"`
    Replicate     int         `schmea:"replicate"`
    EditStop      int         `schema:"edit_stop"`
    JuncEnd       int         `schema:"junc_end"`
    Offset        int         `schema:"offset"`
    Limit         int         `schema:"limit"`
    HasMutation   bool        `schema:"has_mutation"`
    HasAlt        bool        `schema:"has_alt"`
    All           bool        `schema:"all"`
    GrnaEdit      []int       `schema:"grna_edit"`
    GrnaJunc      []int       `schema:"grna_junc"`
}

func (fields *SearchFields) HasKeyMatch(k string) bool {
    key := strings.Split(string(k), ";")
    replicate,_ := strconv.Atoi(key[2])

    if fields.Replicate > 0 && replicate != fields.Replicate {
        return false
    }

    if len(fields.Sample) > 0 {
        match := false
        for _, s := range(fields.Sample) {
            if s == key[1] {
                match = true
                break
            }
        }
        if !match {
            return false
        }
    }
    if fields.Replicate > 0 && fields.Replicate != replicate {
        return false
    }

    return true
}

func (fields *SearchFields) HasMatch(a *Alignment) bool {
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

func NewStorage(dbpath string) (*Storage, error) {
    db, err := bolt.Open(dbpath, 0600, &bolt.Options{Timeout: 1 * time.Second})
    if err != nil {
        return nil, fmt.Errorf("Failed to open database %s - %s", dbpath, err)
    }

    return &Storage{DB: db}, nil
}

func (s *Storage) Search(fields *SearchFields, f func(k string, a *Alignment)) (error) {
    count := 0
    offset := 0

    prefix := ""
    if len(fields.Gene) > 0 {
        prefix +=  fields.Gene
        if len(fields.Sample) == 1 {
            prefix += ";"+fields.Sample[0]
            if fields.Replicate > 0 {
                prefix =  fmt.Sprintf("%s;%d", prefix, fields.Replicate)
            }
        }
    }


    if len(prefix) > 0 {
        err := s.DB.View(func(tx *bolt.Tx) error {
            c := tx.Bucket([]byte("alignments")).Cursor()
            for k, v := c.Seek([]byte(prefix)); bytes.HasPrefix(k, []byte(prefix)); k, v = c.Next() {
                key := string(k)
                if !fields.HasKeyMatch(key) {
                    continue
                }
                a := new(Alignment)
                err := a.UnmarshalBinary(v)
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

            return nil
        })

        return err
    }

    err := s.DB.View(func(tx *bolt.Tx) error {
        b := tx.Bucket([]byte("alignments"))

        b.ForEach(func(k, v []byte) error {
            key := string(k)
            if !fields.HasKeyMatch(key) {
                return nil
            }

            a := new(Alignment)
            err := a.UnmarshalBinary(v)
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

            f(key, a)
            count++
            offset++

            return nil
        })
        return nil
    })

    return err
}

func (s *Storage) PutTemplate(gene string, tmpl *Template) (error) {
    err := s.DB.Update(func(tx *bolt.Tx) error {
        b := tx.Bucket([]byte(BUCKET_TEMPLATES))

        data, err := tmpl.Bytes()
        if err != nil {
            return err
        }

        err = b.Put([]byte(gene), data)
        return err
    })

    return err
}

func (s *Storage) GetTemplate(gene string) (*Template, error) {
    var tmpl *Template
    err := s.DB.View(func(tx *bolt.Tx) error {
        b := tx.Bucket([]byte(BUCKET_TEMPLATES))
        v := b.Get([]byte(gene))

        t, err := NewTemplateFromBytes(v)
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
