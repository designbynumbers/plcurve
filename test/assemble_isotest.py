import glob
import subprocess
import re
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("pdstor", help="pdstor file name",
                    type=str)
args = parser.parse_args()
x = args.pdstor

print "\nRunning carbonize...\n"

try:
    carbonize_results = subprocess.check_output(["carbonizepd",x])
    re.MULTILINE
    xroot = re.search('closing output file (.+)\.\.\.done',carbonize_results);
    if (xroot == None):
        print "couldn't match output file in \n\n" + carbonize_results;
        sys.exit(1)
    print "\t" + x + "\t->\t" + xroot.group(1) 
    carbonize_outfile_name = xroot.group(1) 
except subprocess.CalledProcessError:
    print "carbonize failed on " + x + "\n"
    print carbonize_results
    sys.exit(1)

print "\nExtracting diagram names ..."
names = [ ]
results_list = re.split("\n",carbonize_results)

for l in results_list:
    m = re.search('wrote (pd_create_.+)',l)
    if m:
        names.append(m.group(1))

names = list(set(names))
names.sort()

print "\tfound " + str(len(names)) + " diagram names "

test_name = "./diagram-isotopy-test-tmp.c"
test_file = open(test_name,'w')

print "Writing c test file..."
    
diagram_file = open(carbonize_outfile_name,'r')
test_file.write("/* machine-generated diagram-isomorphism test c code */\n\n")
test_file.write(diagram_file.read())

test_file.write("bool test_diagrams() {\n")
test_file.write("pd_code_t *pd;\n")

for x in names:
    test_file.write("printf(\"testing diagram " + x + "...\");\n")
    test_file.write("pd = " + x + "();\n")
    test_file.write("if (!all_diagram_isotopies_test(pd,\"" + x + "\")) {\n"
                    "  return false;\n"
                    "}\n")
    test_file.write("pd_code_free(&pd);\n\n")

test_file.write("return true;\n}\n\n")
test_file.close()

print "done\n"
