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
	"path/filepath"
	"strings"

	"github.com/Sirupsen/logrus"
	"github.com/ubccr/treat"
)

type LoadOptions struct {
	Gene         string
	Sample       string
	KnockDown    string
	EditBase     string
	TemplatePath string
	FastaPath    string
	Replicate    int
	EditOffset   int
	SkipFrags    bool
	ExcludeSnps  bool
	Force        bool
	Tetracycline bool
	Collapse     bool
}

func cleanName(name string) string {
	// No spaces
	name = strings.Replace(name, " ", "_", -1)
	// No ';'
	name = strings.Replace(name, ";", "_", -1)

	return name
}

func Load(dbpath string, options *LoadOptions) {
	if len(options.Gene) == 0 {
		logrus.Fatal("Gene name is required")
	}
	if len(options.TemplatePath) == 0 {
		logrus.Fatal("Please provide path to templates file")
	}
	if len(options.FastaPath) == 0 {
		logrus.Fatal("Please provide a fasta file to load")
	}
	if len(options.EditBase) != 1 {
		logrus.Fatal("Please provide the edit base")
	}

	if len(options.Sample) == 0 {
		fname := filepath.Base(options.FastaPath)
		options.Sample = fname[:len(fname)-len(filepath.Ext(options.FastaPath))]
	}

	options.Gene = cleanName(options.Gene)
	options.Sample = cleanName(options.Sample)
	options.KnockDown = cleanName(options.KnockDown)

	tmpl, err := treat.NewTemplateFromFasta(options.TemplatePath, treat.FORWARD, rune(options.EditBase[0]))
	if err != nil {
		logrus.Fatalln(err)
	}

	tmpl.SetOffset(options.EditOffset)

	logrus.Printf("Using template Edit Stop Site: %d", tmpl.EditStop)
	logrus.Printf("Using Edit Site numbering offset: %d", tmpl.EditOffset)

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

	_, err = storage.ImportSample(options.FastaPath, options)
	if err != nil {
		logrus.Fatal(err)
	}
}
