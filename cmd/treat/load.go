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
	"bufio"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/Sirupsen/logrus"
	"github.com/aebruno/gofasta"
	"github.com/ubccr/treat"
)

type LoadOptions struct {
	Gene         string
	EditBase     string
	TemplatePath string
	FragmentPath []string
	Primer5      string
	Primer3      string
	Replicate    int
	SkipFrags    bool
	ExcludeSnps  bool
	Force        bool
	Collapse     bool
}


func collapse(tmpl *treat.Template, options *LoadOptions) {
	logrus.Info("Collapsing fragment files excluding primer regions")
	collapsedFiles := make([]string, 0)
	for _, path := range options.FragmentPath {
		logrus.Printf("Processing file: %s", path)
		inFile, err := os.Open(path)
		if err != nil {
			logrus.Fatal(err)
		}
		defer inFile.Close()

		mergeFrags := make(map[string][]*treat.Fragment)
		mergeCounts := make(map[string]uint32)

		for rec := range gofasta.SimpleParser(inFile) {
			frag := treat.NewFragment(rec.Id, rec.Seq, treat.FORWARD, rune(options.EditBase[0]))
			seqNoPrimer, err := frag.StripPrimer(tmpl.Primer5, tmpl.Primer3)
			if err != nil {
				logrus.Fatal(err)
			}
			mergeFrags[seqNoPrimer] = append(mergeFrags[seqNoPrimer], frag)
			mergeCounts[seqNoPrimer] += frag.ReadCount
		}

		outPath := path[:len(path)-len(filepath.Ext(path))]
		outPath += "-primer-collapsed.fasta"
		outFile, err := os.Create(outPath)
		if err != nil {
			logrus.Fatal(err)
		}

		defer outFile.Close()
		writer := bufio.NewWriter(outFile)

		i := 1
		for seq, frags := range mergeFrags {
			max := frags[0].ReadCount
			frag := frags[0]

			for _, f := range frags {
				if f.ReadCount > max {
					max = f.ReadCount
					frag = f
				}
			}

			writer.WriteString(">" + strconv.Itoa(i) + "-" + strconv.Itoa(int(mergeCounts[seq])))
			writer.WriteString("\n")
			writer.WriteString(frag.String())
			writer.WriteString("\n")
			i++
		}
		writer.Flush()
		collapsedFiles = append(collapsedFiles, outPath)
	}

	// Process new collapsed files
	options.FragmentPath = collapsedFiles
}

func Load(dbpath string, options *LoadOptions) {
	if len(options.Gene) == 0 {
		logrus.Fatal("Gene name is required")
	}
	if len(options.TemplatePath) == 0 {
		logrus.Fatal("Please provide path to templates file")
	}
	if options.FragmentPath == nil || len(options.FragmentPath) == 0 {
		logrus.Fatal("Please provide 1 or more fragment files to load")
	}
	if len(options.EditBase) != 1 {
		logrus.Fatal("Please provide the edit base")
	}

    // Clean up gene names
	options.Gene = strings.Replace(options.Gene, " ", "_", -1)

	tmpl, err := treat.NewTemplateFromFasta(options.TemplatePath, treat.FORWARD, rune(options.EditBase[0]))
	if err != nil {
		logrus.Fatalln(err)
	}

    err = tmpl.SetPrimers(options.Primer5, options.Primer3)
    if err != nil {
        logrus.Fatalln(err)
    }

	if options.Collapse {
		collapse(tmpl, options)
	}

	logrus.Printf("Using template Edit Stop Site: %d", tmpl.EditStop)

    storage, err := NewStorageWrite(dbpath)
	if err != nil {
		logrus.Fatal(err)
	}

    err = storage.Initialize()
	if err != nil {
		logrus.Fatal(err)
	}

    err = storage.PutTemplate(options.Gene, tmpl)
	if err != nil {
		logrus.Fatal(err)
	}

	for _, path := range options.FragmentPath {
        _, err := storage.ImportSample(path, options)
        if err != nil {
            logrus.Fatal(err)
        }
	}
}
