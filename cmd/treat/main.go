// Copyright 2015 TREAT Authors. All rights reserved.
// Use of this source code is governed by a BSD style
// license that can be found in the LICENSE file.

package main

import (
    "os"
    "log"
    "path/filepath"
    "github.com/codegangsta/cli"
)

func main() {
    app := cli.NewApp()
    app.Name    = "treat"
    app.Authors = []cli.Author{cli.Author{Name: "Andrew E. Bruno", Email: "aebruno2@buffalo.edu"}}
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
                &cli.StringFlag{Name: "template, t", Usage: "Path to templates file in FASTA format"},
                &cli.StringSliceFlag{Name: "fragment, f", Value: &cli.StringSlice{}, Usage: "One or more fragment FASTA files"},
                &cli.StringFlag{Name: "primer5", Usage: "5' primer sequence"},
                &cli.StringFlag{Name: "primer3", Usage: "3' primer sequence"},
                &cli.StringFlag{Name: "base, b", Value: "T", Usage: "Edit base"},
                &cli.StringFlag{Name: "dir", Usage: "Directory of fragment FASTA files"},
                &cli.Float64Flag{Name: "normalize, n", Value: float64(0), Usage: "Normalize to read count"},
                &cli.BoolFlag{Name: "skip-fragments", Usage: "Do not store raw fragments. Only alignment summary data."},
                &cli.BoolFlag{Name: "exclude-snps", Usage: "Exclude fragments containing SNPs."},
                &cli.BoolFlag{Name: "fastx", Usage: "Parse FASTX/Collpaser header lines"},
                &cli.BoolFlag{Name: "force", Usage: "Force delete gene data if already exists"},
                &cli.BoolFlag{Name: "collapse", Usage: "Collapse fragments excluding primer regions  into a single sequence (while maintaining reads counts)"},
            },
            Action: func(c *cli.Context) {
                dir := c.String("dir")
                fragPaths := c.StringSlice("fragment")
                if fragPaths != nil && len(fragPaths) > 0 && len(dir) > 0 {
                    log.Fatalln("ERROR Use --fragment or --dir not both.")
                }

                if len(dir) > 0 {
                    fa, err := filepath.Glob(dir + "/*.fasta")
                    fasta, err := filepath.Glob(dir + "/*.fa")
                    if err != nil {
                        log.Fatalf("ERROR Failed to find *.fa(sta) files in dir: %s", dir)
                    }
                    fgmts := append(fa, fasta...)

                    for _, f := range fgmts {
                        abs, _ := filepath.Abs(f)
                        fragPaths = append(fragPaths, abs)
                    }
                }

                Load(c.GlobalString("db"), &LoadOptions{
                    Gene:         c.String("gene"),
                    TemplatePath: c.String("template"),
                    FragmentPath: fragPaths,
                    Primer5:      c.String("primer5"),
                    Primer3:      c.String("primer3"),
                    Norm:         c.Float64("normalize"),
                    EditBase:     c.String("base"),
                    SkipFrags:    c.Bool("skip-fragments"),
                    ExcludeSnps:  c.Bool("exclude-snps"),
                    Fastx:        c.Bool("fastx"),
                    Force:        c.Bool("force"),
                    Collapse:        c.Bool("collapse"),
                })
            },
        },
        {
            Name: "align",
            Usage: "Align one or more fragments",
            Flags: []cli.Flag{
                &cli.StringFlag{Name: "template, t", Usage: "Path to templates file in FASTA format"},
                &cli.StringFlag{Name: "fragment, f", Usage: "Path to fragment FASTA file"},
                &cli.StringFlag{Name: "base, b", Value: "T", Usage: "Edit base"},
                &cli.StringFlag{Name: "s1, 1", Usage: "first sequence to align"},
                &cli.StringFlag{Name: "s2, 2", Usage: "second sequence to align"},
                &cli.StringFlag{Name: "primer5", Usage: "5' primer sequence"},
                &cli.StringFlag{Name: "primer3", Usage: "3' primer sequence"},
            },
            Action: func(c *cli.Context) {
                Align(&AlignOptions{
                    TemplatePath: c.String("template"),
                    FragmentPath: c.String("fragment"),
                    Primer5:      c.String("primer5"),
                    Primer3:      c.String("primer3"),
                    EditBase:     c.String("base"),
                    S1:           c.String("s1"),
                    S2:           c.String("s2"),
                })
            },
        },
        {
            Name: "mutant",
            Usage: "Indel mutation analysis",
            Flags: []cli.Flag{
                &cli.StringFlag{Name: "template, t", Usage: "Path to templates file in FASTA format"},
                &cli.StringSliceFlag{Name: "fragment, f", Value: &cli.StringSlice{}, Usage: "One or more fragment FASTA files"},
                &cli.StringFlag{Name: "base, b", Value: "T", Usage: "Edit base"},
                &cli.StringFlag{Name: "primer5", Usage: "5' primer sequence"},
                &cli.StringFlag{Name: "primer3", Usage: "3' primer sequence"},
                &cli.IntFlag{Name: "n", Value: 5, Usage: "Max number of indels to ouptut"},
            },
            Action: func(c *cli.Context) {
                Mutant(&AlignOptions{
                    TemplatePath: c.String("template"),
                    Primer5:      c.String("primer5"),
                    Primer3:      c.String("primer3"),
                    EditBase:     c.String("base"),
                }, c.StringSlice("fragment"), c.Int("n"))
            },
        },
        {
            Name: "server",
            Usage: "Run http server",
            Flags: []cli.Flag{
                &cli.StringFlag{Name: "templates, t", Usage: "Path to html templates directory"},
                &cli.IntFlag{Name: "port, p", Value: 8080, Usage: "Port to listen on"},
            },
            Action: func(c *cli.Context) {
                Server(c.GlobalString("db"), c.String("templates"), c.Int("port"))
            },
        },
        {
            Name: "stats",
            Usage: "Print database stats",
            Flags: []cli.Flag{
                &cli.StringFlag{Name: "gene, g", Usage: "Filter by gene"},
            },
            Action: func(c *cli.Context) {
                Stats(c.GlobalString("db"), c.String("gene"))
            },
        },
        {
            Name: "search",
            Usage: "Search database",
            Flags: []cli.Flag{
                &cli.StringFlag{Name: "gene, g", Usage: "Gene Name"},
                &cli.StringSliceFlag{Name: "sample, s", Value: &cli.StringSlice{}, Usage: "One or more samples"},
                &cli.IntFlag{Name: "edit-stop", Value: -1, Usage: "Edit stop"},
                &cli.IntFlag{Name: "junc-end", Value: -1, Usage: "Junction end"},
                &cli.IntFlag{Name: "junc-len", Value: -1, Usage: "Junction len"},
                &cli.IntFlag{Name: "alt", Value: 0, Usage: "Alt editing region"},
                &cli.IntFlag{Name: "offset,o", Value: 0, Usage: "offset"},
                &cli.IntFlag{Name: "limit,l", Value: 0, Usage: "limit"},
                &cli.BoolFlag{Name: "has-mutation", Usage: "Has mutation"},
                &cli.BoolFlag{Name: "all,a", Usage: "Include all sequences"},
                &cli.BoolFlag{Name: "has-alt", Usage: "Has Alternative Editing"},
                &cli.BoolFlag{Name: "csv", Usage: "Output in csv format"},
                &cli.BoolFlag{Name: "no-header, x", Usage: "Exclude header from output"},
            },
            Action: func(c *cli.Context) {
                Search(c.GlobalString("db"), &SearchFields{
                    Gene:         c.String("gene"),
                    Sample:       c.StringSlice("sample"),
                    EditStop:     c.Int("edit-stop"),
                    JuncLen:      c.Int("junc-len"),
                    JuncEnd:      c.Int("junc-end"),
                    Offset:       c.Int("offset"),
                    Limit:        c.Int("limit"),
                    AltRegion:    c.Int("alt"),
                    HasMutation:  c.Bool("has-mutation"),
                    HasAlt:       c.Bool("has-alt"),
                    All:          c.Bool("all"),
                }, c.Bool("csv"), c.Bool("no-header"))
            },
        }}

    app.Run(os.Args)
}
