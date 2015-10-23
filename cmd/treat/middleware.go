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
	"net/http"

	"github.com/Sirupsen/logrus"
	"github.com/gorilla/context"
)

// DbContext sets the database context for the request
func DbContext(app *Application) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		session, _ := app.cookieStore.Get(r, TREAT_COOKIE_SESSION)
		dbname := session.Values[TREAT_COOKIE_DB]

		if dbname == nil {
			dbname = app.defaultDb
		}

		db, err := app.GetDb(dbname.(string))
		if err != nil {
			logrus.WithFields(logrus.Fields{
				"dbname": dbname,
			}).Warn("Invalid database name. Using default")
			dbname = app.defaultDb
		}

		context.Set(r, "db", db)
	})
}
