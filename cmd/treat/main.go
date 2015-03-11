package main

import (
    "os"
    "github.com/codegangsta/cli"
)

func main() {
    app := cli.NewApp()
    app.Name    = "treat"
    app.Usage   = "Trypanosome RNA Editing Alignment Tool"
    app.Version = "0.0.1"
    app.Flags   = []cli.Flag{
        &cli.StringFlag{Name: "db", Value: "treat.db", Usage: "Path to database file"},
    }
    app.Commands = []cli.Command {
        {
            Name: "load",
            Usage: "Load samples into database",
            Flags: []cli.Flag{
                &cli.StringFlag{Name: "gene, g", Usage: "Gene Name"},
                &cli.StringFlag{Name: "templates, t", Usage: "Path to templates file in FASTA format"},
                &cli.StringSliceFlag{Name: "fragments, f", Value: &cli.StringSlice{}, Usage: "One or more fragment FASTA files"},
                &cli.IntFlag{Name: "replicate, r", Value: 0, Usage: "Biological replicate number"},
                &cli.IntFlag{Name: "primer5", Value: 0, Usage: "5' primer region"},
                &cli.IntFlag{Name: "primer3", Value: 0, Usage: "3' primer region"},
                &cli.Float64Flag{Name: "normalize, n", Usage: "Normalize to read count"},
            },
            Action: func(c *cli.Context) {
                Load(
                    c.GlobalString("db"),
                    c.String("gene"),
                    c.String("templates"),
                    c.StringSlice("fragments"),
                    c.Int("primer5"),
                    c.Int("primer3"),
                    c.Int("replicate"),
                    c.Float64("normalize"))
            },
        },
        {
            Name: "search",
            Usage: "Search database",
            Flags: []cli.Flag{
                &cli.StringFlag{Name: "gene, g", Usage: "Gene Name"},
                &cli.StringFlag{Name: "sample, s", Usage: "Sample Name"},
                &cli.IntFlag{Name: "replicate, r", Value: 0, Usage: "Biological replicate number"},
                &cli.IntFlag{Name: "edit-stop", Value: 0, Usage: "Edit stop"},
                &cli.IntFlag{Name: "junc-end", Value: 0, Usage: "Junction end"},
                &cli.BoolFlag{Name: "has-mutation", Usage: "Has mutation"},
            },
            Action: func(c *cli.Context) {
                Search(c.GlobalString("db"), &SearchFields{
                    Gene:         c.String("gene"),
                    Sample:       c.String("sample"),
                    Replicate:    c.Int("replicate"),
                    EditStop:     c.Int("edit-stop"),
                    JuncEnd:      c.Int("junc-end"),
                    HasMutation:  c.Bool("has-mutation"),
                })
            },
        }}

    app.Run(os.Args)
}
