package treat

import (
    "os"
    "io"
    "fmt"
    "encoding/csv"
    "strconv"
)

type Grna struct {
    Name    string
    Start   uint64
    End     uint64
}

func GrnaFromFile(path string) ([]*Grna, error) {
    f, err := os.Open(path)
    if err != nil {
        return nil, fmt.Errorf("Invalid grna file: %s", err)
    }
    defer f.Close()

    grna := make([]*Grna, 0)

    reader := csv.NewReader(f)

    for {
        record, err := reader.Read()

        if err == io.EOF {
            break
        } else if err != nil {
            return nil, fmt.Errorf("Error parsing file: %s", err)
        }

        if len(record) != 3 {
            return nil, fmt.Errorf("Invalid grna record len. Must be exactly 3 columns but got: %d", len(record))
        }

        start, err := strconv.Atoi(record[1])
        if err != nil {
            return nil, fmt.Errorf("Error start of grna: %s %s", record[0], record[1])
        }

        end, err := strconv.Atoi(record[2])
        if err != nil {
            return nil, fmt.Errorf("Error end of grna: %s %s", record[0], record[2])
        }

        grna = append(grna, &Grna{Name: record[0], Start: uint64(start), End: uint64(end)})
    }

    return grna, nil
}
