package main

import (
    "os"
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
                &cli.StringFlag{Name: "grna", Usage: "Path to grna file"},
                &cli.StringSliceFlag{Name: "fragment, f", Value: &cli.StringSlice{}, Usage: "One or more fragment FASTA files"},
                &cli.IntFlag{Name: "primer5", Value: 0, Usage: "5' primer region"},
                &cli.IntFlag{Name: "primer3", Value: 0, Usage: "3' primer region"},
                &cli.StringFlag{Name: "base, b", Value: "T", Usage: "Edit base"},
                &cli.Float64Flag{Name: "normalize, n", Value: float64(0), Usage: "Normalize to read count"},
                &cli.BoolFlag{Name: "skip-fragments", Usage: "Do not store raw fragments. Only alignment summary data."},
            },
            Action: func(c *cli.Context) {
                Load(c.GlobalString("db"), &LoadOptions{
                    Gene:         c.String("gene"),
                    TemplatePath: c.String("template"),
                    FragmentPath: c.StringSlice("fragment"),
                    Primer5:      c.Int("primer5"),
                    Primer3:      c.Int("primer3"),
                    Norm:         c.Float64("normalize"),
                    GrnaPath:     c.String("grna"),
                    EditBase:     c.String("base"),
                    SkipFrags:    c.Bool("skip-fragments"),
                })
            },
        },
        {
            Name: "align",
            Usage: "Align one or more fragments",
            Flags: []cli.Flag{
                &cli.StringFlag{Name: "template, t", Usage: "Path to templates file in FASTA format"},
                &cli.StringFlag{Name: "grna", Usage: "Path to grna file"},
                &cli.StringFlag{Name: "fragment, f", Usage: "Path to fragment FASTA file"},
                &cli.StringFlag{Name: "base, b", Value: "T", Usage: "Edit base"},
                &cli.StringFlag{Name: "s1, 1", Usage: "first sequence to align"},
                &cli.StringFlag{Name: "s2, 2", Usage: "second sequence to align"},
                &cli.IntFlag{Name: "primer5", Value: 0, Usage: "5' primer region"},
                &cli.IntFlag{Name: "primer3", Value: 0, Usage: "3' primer region"},
            },
            Action: func(c *cli.Context) {
                Align(&AlignOptions{
                    TemplatePath: c.String("template"),
                    FragmentPath: c.String("fragment"),
                    Primer5:      c.Int("primer5"),
                    Primer3:      c.Int("primer3"),
                    GrnaPath:     c.String("grna"),
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
                &cli.IntFlag{Name: "primer5", Value: 0, Usage: "5' primer region"},
                &cli.IntFlag{Name: "primer3", Value: 0, Usage: "3' primer region"},
                &cli.IntFlag{Name: "n", Value: 5, Usage: "Max number of indels to ouptut"},
            },
            Action: func(c *cli.Context) {
                Mutant(&AlignOptions{
                    TemplatePath: c.String("template"),
                    Primer5:      c.Int("primer5"),
                    Primer3:      c.Int("primer3"),
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
                &cli.IntSliceFlag{Name: "grna-edit", Value: &cli.IntSlice{}, Usage: "gRNA over edit stop"},
                &cli.IntSliceFlag{Name: "grna-junc", Value: &cli.IntSlice{}, Usage: "gRNA over junc region"},
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
                    GrnaEdit:     c.IntSlice("grna-edit"),
                    GrnaJunc:     c.IntSlice("grna-junc"),
                }, c.Bool("csv"), c.Bool("no-header"))
            },
        }}

    app.Run(os.Args)
}
