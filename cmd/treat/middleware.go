// Copyright 2015 TREAT Authors. All rights reserved.
// Use of this source code is governed by a BSD style
// license that can be found in the LICENSE file.

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
