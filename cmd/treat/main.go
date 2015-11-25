// Copyright 2015 TREAT Authors. All rights reserved.
//
// This file is part of TREAT.
//
// TREAT is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// TREAT is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with TREAT.  If not, see <http://www.gnu.org/licenses/>.

package main

import (
	"os"
	"path/filepath"

	"github.com/Sirupsen/logrus"
	"github.com/codegangsta/cli"
)

func main() {
	app := cli.NewApp()
	app.Name = "treat"
	app.Copyright = `Copyright 2015 TREAT Authors.  

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.`
	app.Authors = []cli.Author{
		{Name: "Andrew E. Bruno", Email: "aebruno2@buffalo.edu"},
		{Name: "Rachel Simpson", Email: "rachel.simpson64@gmail.com"},
		{Name: "Laurie Read", Email: "lread@buffalo.edu"}}
	app.Usage = "Trypanosome RNA Editing Alignment Tool"
	app.Version = "0.0.1"
	app.Flags = []cli.Flag{
		&cli.StringFlag{Name: "db", Value: "treat.db", Usage: "Path to database file"},
	}
	app.Commands = []cli.Command{
		{
			Name:  "load",
			Usage: "Load samples into database",
			Flags: []cli.Flag{
				&cli.StringFlag{Name: "gene, g", Usage: "Gene Name"},
				&cli.StringFlag{Name: "template, t", Usage: "Path to templates file in FASTA format"},
				&cli.StringSliceFlag{Name: "fragment, f", Value: &cli.StringSlice{}, Usage: "One or more fragment FASTA files"},
				&cli.StringFlag{Name: "primer5", Usage: "5' primer sequence"},
				&cli.StringFlag{Name: "primer3", Usage: "3' primer sequence"},
				&cli.StringFlag{Name: "base, b", Value: "T", Usage: "Edit base"},
				&cli.StringFlag{Name: "dir", Usage: "Directory of fragment FASTA files"},
				&cli.BoolFlag{Name: "skip-fragments", Usage: "Do not store raw fragments. Only alignment summary data."},
				&cli.BoolFlag{Name: "exclude-snps", Usage: "Exclude fragments containing SNPs."},
				&cli.BoolFlag{Name: "force", Usage: "Force delete gene data if already exists"},
				&cli.BoolFlag{Name: "collapse", Usage: "Collapse fragments excluding primer regions  into a single sequence (while maintaining reads counts)"},
			},
			Action: func(c *cli.Context) {
				dir := c.String("dir")
				fragPaths := c.StringSlice("fragment")
				if fragPaths != nil && len(fragPaths) > 0 && len(dir) > 0 {
					logrus.Fatal("Use --fragment or --dir not both.")
				}

				if len(dir) > 0 {
					fa, _ := filepath.Glob(filepath.Join(dir, "*.fasta"))
					fasta, _ := filepath.Glob(filepath.Join(dir, "*.fa"))
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
					EditBase:     c.String("base"),
					SkipFrags:    c.Bool("skip-fragments"),
					ExcludeSnps:  c.Bool("exclude-snps"),
					Force:        c.Bool("force"),
					Collapse:     c.Bool("collapse"),
				})
			},
		},
		{
			Name:  "align",
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
			Name:  "mutant",
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
			Name:  "server",
			Usage: "Run http server",
			Flags: []cli.Flag{
				&cli.StringFlag{Name: "templates, t", Usage: "Path to html templates directory"},
				&cli.IntFlag{Name: "port, p", Value: 8080, Usage: "Port to listen on"},
				&cli.BoolFlag{Name: "enable-cache", Usage: "Enable url caching"},
			},
			Action: func(c *cli.Context) {
				Server(c.GlobalString("db"), c.String("templates"), c.Int("port"), c.Bool("enable-cache"))
			},
		},
		{
			Name:  "stats",
			Usage: "Print database stats",
			Flags: []cli.Flag{
				&cli.StringFlag{Name: "gene, g", Usage: "Filter by gene"},
				&cli.BoolFlag{Name: "unique, u", Usage: "Use unique fragment counts only"},
				&cli.BoolFlag{Name: "norm, n", Usage: "Use normalized fragment counts only"},
			},
			Action: func(c *cli.Context) {
				ShowStats(c.GlobalString("db"), c.String("gene"), c.Bool("unique"), c.Bool("norm"))
			},
		},
		{
			Name:  "norm",
			Usage: "Normalize read counts",
			Flags: []cli.Flag{
				&cli.StringFlag{Name: "gene, g", Usage: "Gene name (all by default)"},
				&cli.Float64Flag{Name: "normalize, n", Value: float64(0), Usage: "Normalize to read count"},
			},
			Action: func(c *cli.Context) {
				Normalize(c.GlobalString("db"), c.String("gene"), c.Float64("normalize"))
			},
		},
		{
			Name:  "search",
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
					Gene:        c.String("gene"),
					Sample:      c.StringSlice("sample"),
					EditStop:    c.Int("edit-stop"),
					JuncLen:     c.Int("junc-len"),
					JuncEnd:     c.Int("junc-end"),
					Offset:      c.Int("offset"),
					Limit:       c.Int("limit"),
					AltRegion:   c.Int("alt"),
					HasMutation: c.Bool("has-mutation"),
					HasAlt:      c.Bool("has-alt"),
					All:         c.Bool("all"),
				}, c.Bool("csv"), c.Bool("no-header"))
			},
		}}

	app.Run(os.Args)
}
