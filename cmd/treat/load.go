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
	"strings"

	"github.com/Sirupsen/logrus"
	"github.com/ubccr/treat"
)

type LoadOptions struct {
	Gene         string
	EditBase     string
	TemplatePath string
	FragmentPath []string
	Replicate    int
	SkipFrags    bool
	ExcludeSnps  bool
	Force        bool
	Collapse     bool
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
