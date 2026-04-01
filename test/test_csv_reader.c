/*
 * test_csv_reader.c -- Tests for plc_read_csv and plc_read_tsv
 *
 * Copyright 2024 The University of Georgia
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "plCurve.h"

static int failures = 0;

#define ASSERT_TRUE(cond, msg) do { \
  if (!(cond)) { \
    fprintf(stderr, "FAIL: %s (line %d)\n", msg, __LINE__); \
    failures++; \
  } \
} while (0)

#define ASSERT_NEAR(a, b, tol, msg) do { \
  if (fabs((a) - (b)) > (tol)) { \
    fprintf(stderr, "FAIL: %s: expected %g, got %g (line %d)\n", msg, (double)(b), (double)(a), __LINE__); \
    failures++; \
  } \
} while (0)

/* Helper: write a string to a temp file and return it rewound */
static FILE *string_to_file(const char *content)
{
  FILE *f = tmpfile();
  if (f == NULL) return NULL;
  fputs(content, f);
  rewind(f);
  return f;
}

static void test_csv_closed_triangle(void)
{
  const char *data = "0.0, 0.0, 0.0\n1.0, 0.0, 0.0\n0.5, 1.0, 0.0\n";
  FILE *f = string_to_file(data);
  int err = 0;
  char errmsg[256] = "";

  plCurve *L = plc_read_csv(f, false, &err, errmsg, sizeof(errmsg));
  fclose(f);

  ASSERT_TRUE(L != NULL, "csv closed triangle: L not null");
  ASSERT_TRUE(err == 0, "csv closed triangle: no error");
  if (L != NULL) {
    ASSERT_TRUE(L->nc == 1, "csv closed triangle: 1 component");
    ASSERT_TRUE(L->cp[0].nv == 3, "csv closed triangle: 3 vertices");
    ASSERT_TRUE(L->cp[0].open == false, "csv closed triangle: closed");
    ASSERT_NEAR(L->cp[0].vt[0].c[0], 0.0, 1e-12, "csv closed triangle: v0.x");
    ASSERT_NEAR(L->cp[0].vt[1].c[0], 1.0, 1e-12, "csv closed triangle: v1.x");
    ASSERT_NEAR(L->cp[0].vt[2].c[1], 1.0, 1e-12, "csv closed triangle: v2.y");
    plc_free(L);
  }
}

static void test_tsv_closed_triangle(void)
{
  const char *data = "0.0\t0.0\t0.0\n1.0\t0.0\t0.0\n0.5\t1.0\t0.0\n";
  FILE *f = string_to_file(data);
  int err = 0;
  char errmsg[256] = "";

  plCurve *L = plc_read_tsv(f, false, &err, errmsg, sizeof(errmsg));
  fclose(f);

  ASSERT_TRUE(L != NULL, "tsv closed triangle: L not null");
  ASSERT_TRUE(err == 0, "tsv closed triangle: no error");
  if (L != NULL) {
    ASSERT_TRUE(L->nc == 1, "tsv closed triangle: 1 component");
    ASSERT_TRUE(L->cp[0].nv == 3, "tsv closed triangle: 3 vertices");
    ASSERT_TRUE(L->cp[0].open == false, "tsv closed triangle: closed");
    ASSERT_NEAR(L->cp[0].vt[1].c[0], 1.0, 1e-12, "tsv closed triangle: v1.x");
    plc_free(L);
  }
}

static void test_csv_open_curve(void)
{
  const char *data = "0.0,0.0,0.0\n1.0,0.0,0.0\n2.0,0.0,0.0\n";
  FILE *f = string_to_file(data);
  int err = 0;
  char errmsg[256] = "";

  plCurve *L = plc_read_csv(f, true, &err, errmsg, sizeof(errmsg));
  fclose(f);

  ASSERT_TRUE(L != NULL, "csv open: L not null");
  ASSERT_TRUE(err == 0, "csv open: no error");
  if (L != NULL) {
    ASSERT_TRUE(L->cp[0].open == true, "csv open: is open");
    ASSERT_TRUE(L->cp[0].nv == 3, "csv open: 3 vertices");
    plc_free(L);
  }
}

static void test_csv_duplicate_removed_closed(void)
{
  /* Last vertex duplicates the first -- should be removed for closed curve */
  const char *data = "0.0,0.0,0.0\n1.0,0.0,0.0\n0.5,1.0,0.0\n0.0,0.0,0.0\n";
  FILE *f = string_to_file(data);
  int err = 0;
  char errmsg[256] = "";

  plCurve *L = plc_read_csv(f, false, &err, errmsg, sizeof(errmsg));
  fclose(f);

  ASSERT_TRUE(L != NULL, "csv dup closed: L not null");
  ASSERT_TRUE(err == 0, "csv dup closed: no error");
  if (L != NULL) {
    ASSERT_TRUE(L->cp[0].nv == 3, "csv dup closed: 3 verts (dup removed)");
    plc_free(L);
  }
}

static void test_csv_duplicate_kept_open(void)
{
  /* Last vertex duplicates the first -- should NOT be removed for open curve */
  const char *data = "0.0,0.0,0.0\n1.0,0.0,0.0\n0.5,1.0,0.0\n0.0,0.0,0.0\n";
  FILE *f = string_to_file(data);
  int err = 0;
  char errmsg[256] = "";

  plCurve *L = plc_read_csv(f, true, &err, errmsg, sizeof(errmsg));
  fclose(f);

  ASSERT_TRUE(L != NULL, "csv dup open: L not null");
  ASSERT_TRUE(err == 0, "csv dup open: no error");
  if (L != NULL) {
    ASSERT_TRUE(L->cp[0].nv == 4, "csv dup open: 4 verts (dup kept)");
    plc_free(L);
  }
}

static void test_csv_non_duplicate_kept(void)
{
  /* Last vertex does not duplicate the first -- should be kept */
  const char *data = "0.0,0.0,0.0\n1.0,0.0,0.0\n0.5,1.0,0.0\n0.0,0.0,1.0\n";
  FILE *f = string_to_file(data);
  int err = 0;
  char errmsg[256] = "";

  plCurve *L = plc_read_csv(f, false, &err, errmsg, sizeof(errmsg));
  fclose(f);

  ASSERT_TRUE(L != NULL, "csv non-dup: L not null");
  ASSERT_TRUE(err == 0, "csv non-dup: no error");
  if (L != NULL) {
    ASSERT_TRUE(L->cp[0].nv == 4, "csv non-dup: 4 verts kept");
    plc_free(L);
  }
}

static void test_csv_spaces_around_delimiters(void)
{
  /* Mimics randompolygon output: "%12g, %12g, %12g \n" */
  const char *data = "     0.5 ,      1.0 ,      0.0 \n     1.0 ,      0.0 ,      0.0 \n     0.0 ,      0.0 ,      1.0 \n";
  FILE *f = string_to_file(data);
  int err = 0;
  char errmsg[256] = "";

  plCurve *L = plc_read_csv(f, false, &err, errmsg, sizeof(errmsg));
  fclose(f);

  ASSERT_TRUE(L != NULL, "csv spaces: L not null");
  ASSERT_TRUE(err == 0, "csv spaces: no error");
  if (L != NULL) {
    ASSERT_TRUE(L->cp[0].nv == 3, "csv spaces: 3 vertices");
    ASSERT_NEAR(L->cp[0].vt[0].c[0], 0.5, 1e-12, "csv spaces: v0.x");
    ASSERT_NEAR(L->cp[0].vt[0].c[1], 1.0, 1e-12, "csv spaces: v0.y");
    plc_free(L);
  }
}

static void test_csv_comments_and_blanks(void)
{
  const char *data =
    "# This is a header\n"
    "\n"
    "0.0, 0.0, 0.0\n"
    "# Another comment\n"
    "1.0, 0.0, 0.0\n"
    "\n"
    "0.5, 1.0, 0.0\n";
  FILE *f = string_to_file(data);
  int err = 0;
  char errmsg[256] = "";

  plCurve *L = plc_read_csv(f, false, &err, errmsg, sizeof(errmsg));
  fclose(f);

  ASSERT_TRUE(L != NULL, "csv comments: L not null");
  ASSERT_TRUE(err == 0, "csv comments: no error");
  if (L != NULL) {
    ASSERT_TRUE(L->cp[0].nv == 3, "csv comments: 3 vertices");
    plc_free(L);
  }
}

static void test_csv_too_few_verts_closed(void)
{
  /* Only 2 vertices -- closed curves need at least 3 */
  const char *data = "0.0,0.0,0.0\n1.0,0.0,0.0\n";
  FILE *f = string_to_file(data);
  int err = 0;
  char errmsg[256] = "";

  plCurve *L = plc_read_csv(f, false, &err, errmsg, sizeof(errmsg));
  fclose(f);

  ASSERT_TRUE(L == NULL, "csv too few closed: L is null");
  ASSERT_TRUE(err == PLC_E_TOO_FEW_VERTS, "csv too few closed: correct error code");
}

static void test_csv_too_few_verts_open(void)
{
  /* Only 1 vertex -- open curves need at least 2 */
  const char *data = "0.0,0.0,0.0\n";
  FILE *f = string_to_file(data);
  int err = 0;
  char errmsg[256] = "";

  plCurve *L = plc_read_csv(f, true, &err, errmsg, sizeof(errmsg));
  fclose(f);

  ASSERT_TRUE(L == NULL, "csv too few open: L is null");
  ASSERT_TRUE(err == PLC_E_TOO_FEW_VERTS, "csv too few open: correct error code");
}

static void test_csv_malformed_line(void)
{
  const char *data = "0.0,0.0,0.0\nhello,world,foo\n0.5,1.0,0.0\n";
  FILE *f = string_to_file(data);
  int err = 0;
  char errmsg[256] = "";

  plCurve *L = plc_read_csv(f, false, &err, errmsg, sizeof(errmsg));
  fclose(f);

  ASSERT_TRUE(L == NULL, "csv malformed: L is null");
  ASSERT_TRUE(err == PLC_E_BAD_CSV_LINE, "csv malformed: correct error code");
}

static void test_csv_open_two_verts(void)
{
  /* Exactly 2 vertices -- should work for open curves */
  const char *data = "0.0,0.0,0.0\n1.0,0.0,0.0\n";
  FILE *f = string_to_file(data);
  int err = 0;
  char errmsg[256] = "";

  plCurve *L = plc_read_csv(f, true, &err, errmsg, sizeof(errmsg));
  fclose(f);

  ASSERT_TRUE(L != NULL, "csv open 2 verts: L not null");
  ASSERT_TRUE(err == 0, "csv open 2 verts: no error");
  if (L != NULL) {
    ASSERT_TRUE(L->cp[0].nv == 2, "csv open 2 verts: 2 vertices");
    ASSERT_TRUE(L->cp[0].open == true, "csv open 2 verts: is open");
    plc_free(L);
  }
}

static void test_csv_dup_removal_drops_below_min(void)
{
  /* 3 vertices but last duplicates first -> 2 verts after removal -> too few for closed */
  const char *data = "0.0,0.0,0.0\n1.0,0.0,0.0\n0.0,0.0,0.0\n";
  FILE *f = string_to_file(data);
  int err = 0;
  char errmsg[256] = "";

  plCurve *L = plc_read_csv(f, false, &err, errmsg, sizeof(errmsg));
  fclose(f);

  ASSERT_TRUE(L == NULL, "csv dup below min: L is null");
  ASSERT_TRUE(err == PLC_E_TOO_FEW_VERTS, "csv dup below min: correct error code");
}

int main(void)
{
  printf("test_csv_reader: running tests...\n");

  test_csv_closed_triangle();
  test_tsv_closed_triangle();
  test_csv_open_curve();
  test_csv_duplicate_removed_closed();
  test_csv_duplicate_kept_open();
  test_csv_non_duplicate_kept();
  test_csv_spaces_around_delimiters();
  test_csv_comments_and_blanks();
  test_csv_too_few_verts_closed();
  test_csv_too_few_verts_open();
  test_csv_malformed_line();
  test_csv_open_two_verts();
  test_csv_dup_removal_drops_below_min();

  if (failures == 0) {
    printf("test_csv_reader: all tests passed.\n");
    return 0;
  } else {
    printf("test_csv_reader: %d test(s) FAILED.\n", failures);
    return 1;
  }
}
