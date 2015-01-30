from collections import defaultdict
def kt_indexer(kt):
    NX_MUL = 10**9
    FACTOR_MUL = 10**4
    NF_MUL = 10**8
    IDX_MUL = 10**2
    MIR_MUL = 10

    val = 0
    for i,factor in enumerate(kt.split("#")):
        val += NF_MUL
        factor = factor.strip()
        if "*" in factor:
            val += MIR_MUL
        factor = factor.split("*")[0]
        #print factor
        nx, idx = factor.split("_")
        val += FACTOR_MUL* int(nx) * (10**(2-i))
        val += NX_MUL * int(nx)
        val += IDX_MUL * int(idx)
    return val

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Join tsv tables into one")
    parser.add_argument('files', metavar="file", type=argparse.FileType('rb'), nargs="+",
                        help="a csv file to join")
    parser.add_argument('output', metavar="output", type=argparse.FileType('wb'),
                        help="a file to create")
    parser.add_argument('--no-symmetry', action='store_true',
                        help="do not use symmetry of data to fill out table")

    args = parser.parse_args()
    table = defaultdict(dict)

    import csv
    for input_csv in args.files:
        input_reader = csv.DictReader(input_csv, delimiter="\t", dialect="excel-tab")
        for row in input_reader:
            for key, val in row.iteritems():
                key = key.replace(" ", "")
                if key and key != "K":
                    #print key, row["K"]
                    print kt_indexer(key), kt_indexer(row["K"])
                    row_key = row["K"].replace(" ", "")
                    #print key, row_key, table[key]
                    assert row_key not in table[key]
                    val = val.replace("-", "~")
                    table[key][row_key] = val
    #print table
    fieldnames = ["K"] + sorted(table.keys(), key=kt_indexer)
    print fieldnames

    table_transpose = defaultdict(dict)
    for rowname, row in table.iteritems():
        for colname, col in row.iteritems():
            table_transpose[colname][rowname] = col

    print table_transpose
    #assert(False)

    rownames = sorted(table_transpose.keys(), key=kt_indexer)
    for rowname in rownames:
        for colname in fieldnames[1:]:
            if colname in table_transpose[rowname] and rowname in table_transpose[colname]:
                assert table_transpose[rowname][colname] == table_transpose[colname][rowname]
            elif rowname in table_transpose[colname]:
                table_transpose[rowname][colname] = table_transpose[colname][rowname]

    out_writer = csv.DictWriter(args.output, fieldnames)
    out_writer.writeheader()
    for rowname in rownames:
        row = table_transpose[rowname]
        outrow = dict(row)
        outrow["K"] = rowname
        out_writer.writerow(outrow)
