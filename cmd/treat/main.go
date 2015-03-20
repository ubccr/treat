package main

import (
    "os"
    "log"
    "github.com/codegangsta/cli"
    "github.com/ubccr/treat"
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
                &cli.StringFlag{Name: "template, t", Usage: "Path to templates file in FASTA format"},
                &cli.StringFlag{Name: "grna", Usage: "Path to grna file"},
                &cli.StringSliceFlag{Name: "fragment, f", Value: &cli.StringSlice{}, Usage: "One or more fragment FASTA files"},
                &cli.IntFlag{Name: "primer5", Value: 0, Usage: "5' primer region"},
                &cli.IntFlag{Name: "primer3", Value: 0, Usage: "3' primer region"},
                &cli.Float64Flag{Name: "normalize, n", Value: float64(0), Usage: "Normalize to read count"},
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
                })
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
                //decoder.IgnoreUnknownKeys(true)
                Server(c.GlobalString("db"), c.String("templates"), c.Int("port"))
            },
        },
        {
            Name: "stats",
            Usage: "Print database stats",
            Action: func(c *cli.Context) {
                s, err := treat.NewStorage(c.GlobalString("db"))
                if err != nil {
                    log.Fatal(err)
                }
                s.Stats()
            },
        },
        {
            Name: "search",
            Usage: "Search database",
            Flags: []cli.Flag{
                &cli.StringFlag{Name: "gene, g", Usage: "Gene Name"},
                &cli.StringSliceFlag{Name: "sample, s", Value: &cli.StringSlice{}, Usage: "One or more samples"},
                &cli.IntFlag{Name: "edit-stop", Value: 0, Usage: "Edit stop"},
                &cli.IntFlag{Name: "junc-end", Value: 0, Usage: "Junction end"},
                &cli.IntFlag{Name: "junc-len", Value: 0, Usage: "Junction len"},
                &cli.IntFlag{Name: "offset,o", Value: 0, Usage: "offset"},
                &cli.IntFlag{Name: "limit,l", Value: 0, Usage: "limit"},
                &cli.BoolFlag{Name: "has-mutation", Usage: "Has mutation"},
                &cli.BoolFlag{Name: "all,a", Usage: "Include all sequences"},
                &cli.BoolFlag{Name: "has-alt", Usage: "Has Alternative Editing"},
                &cli.IntSliceFlag{Name: "grna-edit", Value: &cli.IntSlice{}, Usage: "gRNA over edit stop"},
                &cli.IntSliceFlag{Name: "grna-junc", Value: &cli.IntSlice{}, Usage: "gRNA over junc region"},
            },
            Action: func(c *cli.Context) {
                Search(c.GlobalString("db"), &treat.SearchFields{
                    Gene:         c.String("gene"),
                    Sample:       c.StringSlice("sample"),
                    EditStop:     c.Int("edit-stop"),
                    JuncLen:      c.Int("junc-len"),
                    JuncEnd:      c.Int("junc-end"),
                    Offset:       c.Int("offset"),
                    Limit:        c.Int("limit"),
                    HasMutation:  c.Bool("has-mutation"),
                    HasAlt:       c.Bool("has-alt"),
                    All:          c.Bool("all"),
                    GrnaEdit:     c.IntSlice("grna-edit"),
                    GrnaJunc:     c.IntSlice("grna-junc"),
                })
            },
        }}

    app.Run(os.Args)
}
