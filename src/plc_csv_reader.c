/*
 * plc_csv_reader.c -- Read plCurve data from CSV and TSV files
 *
 * Copyright 2024 The University of Georgia
 *
 * This file is part of plCurve.
 *
 * plCurve is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "plCurve.h"

#define CSV_INITIAL_CAP 256
#define CSV_LINE_MAX 4096

static plCurve *plc_read_delimited(FILE *file, char delim, bool open,
                                   int *error_num, char error_str[],
                                   size_t error_str_len)
{
  char line[CSV_LINE_MAX];
  plc_vector *buf = NULL;
  int nv = 0;
  int cap = CSV_INITIAL_CAP;
  int line_num = 0;
  plCurve *L = NULL;

  *error_num = 0;
  if (error_str != NULL && error_str_len > 0) {
    error_str[0] = '\0';
  }

  buf = (plc_vector *)malloc(cap * sizeof(plc_vector));
  if (buf == NULL) {
    *error_num = PLC_E_BAD_CSV_LINE;
    if (error_str != NULL) {
      snprintf(error_str, error_str_len, "plc_read_csv: out of memory");
    }
    return NULL;
  }

  while (fgets(line, sizeof(line), file) != NULL) {
    line_num++;
    char *p = line;

    /* Skip leading whitespace */
    while (*p == ' ' || *p == '\t') p++;

    /* Skip blank lines and comment lines */
    if (*p == '\n' || *p == '\r' || *p == '\0' || *p == '#') {
      continue;
    }

    /* Parse three doubles separated by delim.
       We use strtod which naturally skips leading whitespace,
       so "  1.0 , 2.0 , 3.0" works fine. We just need to
       find and skip the delimiter between values. */
    double coords[3];
    char *endp;
    int i;

    for (i = 0; i < 3; i++) {
      errno = 0;
      coords[i] = strtod(p, &endp);
      if (endp == p || errno != 0) {
        *error_num = PLC_E_BAD_CSV_LINE;
        if (error_str != NULL) {
          snprintf(error_str, error_str_len,
                   "plc_read_csv: parse error on line %d", line_num);
        }
        free(buf);
        return NULL;
      }
      p = endp;

      if (i < 2) {
        /* Skip whitespace (but not the delimiter itself) before the delimiter */
        while (*p != delim && (*p == ' ' || *p == '\t')) p++;
        if (*p == delim) {
          p++;
        } else {
          *error_num = PLC_E_BAD_CSV_LINE;
          if (error_str != NULL) {
            snprintf(error_str, error_str_len,
                     "plc_read_csv: expected delimiter on line %d", line_num);
          }
          free(buf);
          return NULL;
        }
      }
    }

    /* Grow buffer if needed */
    if (nv >= cap) {
      cap *= 2;
      plc_vector *newbuf = (plc_vector *)realloc(buf, cap * sizeof(plc_vector));
      if (newbuf == NULL) {
        *error_num = PLC_E_BAD_CSV_LINE;
        if (error_str != NULL) {
          snprintf(error_str, error_str_len, "plc_read_csv: out of memory");
        }
        free(buf);
        return NULL;
      }
      buf = newbuf;
    }

    buf[nv] = plc_build_vect(coords[0], coords[1], coords[2]);
    nv++;
  }

  /* If closed and last vertex duplicates first, remove duplicate */
  if (!open && nv >= 2 && plc_distance(buf[0], buf[nv - 1]) < 1e-10) {
    nv--;
  }

  /* Check minimum vertex count */
  int min_verts = open ? 2 : 3;
  if (nv < min_verts) {
    *error_num = PLC_E_TOO_FEW_VERTS;
    if (error_str != NULL) {
      snprintf(error_str, error_str_len,
               "plc_read_csv: need at least %d vertices for %s curve, got %d",
               min_verts, open ? "open" : "closed", nv);
    }
    free(buf);
    return NULL;
  }

  /* Build the plCurve */
  int cc = 0;
  L = plc_new(1, &nv, &open, &cc);
  if (L == NULL) {
    *error_num = PLC_E_BAD_CSV_LINE;
    if (error_str != NULL) {
      snprintf(error_str, error_str_len, "plc_read_csv: plc_new failed");
    }
    free(buf);
    return NULL;
  }

  /* Copy vertices */
  int j;
  for (j = 0; j < nv; j++) {
    L->cp[0].vt[j] = buf[j];
  }

  plc_fix_wrap(L);
  free(buf);

  return L;
}

plCurve *plc_read_csv(FILE *file, bool open, int *error_num,
                      char error_str[], size_t error_str_len)
{
  return plc_read_delimited(file, ',', open, error_num, error_str, error_str_len);
}

plCurve *plc_read_tsv(FILE *file, bool open, int *error_num,
                      char error_str[], size_t error_str_len)
{
  return plc_read_delimited(file, '\t', open, error_num, error_str, error_str_len);
}
