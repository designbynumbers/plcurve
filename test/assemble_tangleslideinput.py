import glob
import subprocess
import re
import sys
import os.path

pdstor_files = glob.glob('*.pdstor')
print ("Tangle Slide Input Test Assembly Script\n"
      "==========================\n"
      "Found pdstor files...\n")
print "\n".join("\t"+str(x) for x in pdstor_files)
print "Looking for corresponding .c files and moving to .doc..."
for x in pdstor_files:
    m = re.match('(.+).pdstor',x)
    cname = m.group(1) + ".c"
    docname = m.group(1) + ".doc"
    if (os.path.isfile(cname)):
        # We want to slurp this and remove the \r characters
        cfile_object = open(cname,'r')
        cfile_contents = cfile_object.read()
        cfile_contents = cfile_contents.replace('\r','')
        cfile_object.close()

        docfile_object = open(docname,'w')
        docfile_object.write(cfile_contents)
        docfile_object.close()
        
        print cname + " -> " + docname
    else:
        print "Couldn\'t match "+ x + "to " + cname;
        sys.exit(1)
        
print "\nRunning carbonize...\n"
for x in pdstor_files:
    try:
        carbonize_results = subprocess.check_output(["carbonizepd",x])
        xroot = re.search('(.+)\.pdstor',x);
        print "\t" + x + "\t->\t" + xroot.group(1) + ".c"
    except subprocess.CalledProcessError:
        print "carbonize failed on " + x + "\n"
        print carbonize_results
        sys.exit(1)

print "\nExtracting test names..."
names = [ ]
for x in pdstor_files:
    m = re.match('tangle_slide_input_test(.+).pdstor',x)
    names.append(m.group(1))
    
names = list(set(names))
names.sort()

print "\tfound test names " + ", ".join(names)
print "\nGenerating combined individual test files ... "

cnames = []
test_names = []

for x in names:
    
    fname = "./tangle_slide_input_test{name}.c".format(name=x)
    # the first thing we want to do is to slurp the contents of this file

    try:
        infile = open(fname,'r')
    except EnvironmentError:
        print "couldn't open " + fname + " for reading"
        sys.exit(1)
        
    pdcode_output = infile.read()
    infile.close()

    print "\tSlurped the c form of the pdstor from " + fname

    doc_name = "./tangle_slide_input_test{name}.doc".format(name=x)

    # now we need to slurp the doc file
    try:
        docfile = open(doc_name,'r')
    except EnvironmentError:
        print "couldn't open " + doc_name + " for reading"
        sys.exit(1)

    docstring = docfile.read()
    docfile.close()

    print "\tSlurped the documentation for the test from " + doc_name

    boilerplatefunction = r"""
bool pdint_check_tslide_data_ok_and_find_te_test{NAME}() {{

  {DOCSTRING}

  printf("--------------------------------------------------\n"
    	 "pdint_check_tslide_data_ok_and_find_te test {NAME}\n"
	     "--------------------------------------------------\n");

  printf("creating pd...");
  pd_code_t *pd = pd_create_tangle_slide_input_test{NAME}_0();
  if (!pd_ok(pd)) {{

    printf("fail (doesn't pass pd_ok)\n");
    return false;

  }}
  
  printf("pass (passes pd_ok)\n");

  printf("creating basic tangle information...");
  pd_tangle_t *t = pd_tangle_new(nedges);

  pd_idx_t i;
  for(i=0;i<t->nedges;i++) {{
  
    t->edge[i] = tangle_edges[i];
    t->face[i] = tangle_faces[i];

  }}

  printf("done (tangle has %d edges)\n",nedges);

  printf("running pd_regenerate_tangle...");
  pd_regenerate_tangle(pd,t);

  if (!pd_tangle_ok(pd,t)) {{

    printf("fail (doesn't pass pd_tangle_ok)\n");
    return false;
    
  }}

  printf("pass (didn't crash, passes pd_tangle_ok)\n");
  printf("running pdint_check_tslide_data_ok_and_find_te ");

  pd_idx_t completely_fake_memory_address;
  pd_idx_t *computed_tangle_slide_edges = &completely_fake_memory_address;
  /* We want to initialize this to something that's not NULL, since we need it
     to be NULL if the input data is invalid. */
  
  bool computed_ts_input;

  if (valid_ts_input) {{

    printf("(expect true)...");
    
  }} else {{

    printf("(expect false)...");

  }}

  bool strand_goes_OVER;

  computed_ts_input = pdint_check_tslide_data_ok_and_find_te(pd,t,noverstrand_edges,
                                                             overstrand_edges,
                                                             border_faces,
                                                             &computed_tangle_slide_edges,&strand_goes_OVER);
                                                             
  if (computed_ts_input != valid_ts_input) {{

     printf("FAIL");

     if (computed_ts_input) {{

        printf(" (got true)\n");

     }} else {{

        printf(" (got false)\n");

     }}

     return false;

  }}

  printf("pass");

  if (computed_ts_input) {{

        printf(" (got true)\n");

  }} else {{

        printf(" (got false)\n");

  }}

  if (computed_ts_input) {{ /* If there ARE tangle edges, try to match them... */

    printf("checking address of computed_tangle_slide_edges...");
    if (computed_tangle_slide_edges == &completely_fake_memory_address) {{

      printf("FAIL (not updated in call)\n");
      return false;

    }}

    if (computed_tangle_slide_edges == NULL) {{

      printf("FAIL (set to NULL, even though input is valid)\n");
      return false;

    }}

    printf("pass (set to a new memory address != NULL)\n");

    printf("checking computed tangle_slide_edges against expected...");

    pd_idx_t i;

    for(i=0;i<noverstrand_edges-1;i++) {{

     if (computed_tangle_slide_edges[i] != tangle_slide_edges[i]) {{

        printf("FAIL.\nExpected and computed tangle edges don't match at pos %d\n",i);

        pd_idx_t j;

        printf("Expected tangle edges: ");
        for(j=0;j<noverstrand_edges;j++) {{ printf("%4d ",tangle_slide_edges[j]); }}
        printf("\nComputed tangle edges: ");
        for(j=0;j<noverstrand_edges;j++) {{ printf("%4d ",computed_tangle_slide_edges[j]); }}
        printf("\n                       ");
        for(j=0;j<i;j++) {{ printf("     "); }}
        printf("-----\n");

        return false;

      }}

    }}

    printf("pass (lists of edges match)\n");

  }} else {{

    printf("checking computed_tangle_slide_edges == NULL ...");

    if (computed_tangle_slide_edges != NULL) {{

       printf("FAIL.\n");
       return false;

    }}

    printf("pass\n");

  }}
  
  printf("housecleaning...");
  
  pd_code_free(&pd);
  pd_tangle_free(&t);
  free(computed_tangle_slide_edges);

  printf("done\n");

  printf("--------------------------------------------------------\n"
         "pdint_check_tslide_data_ok_and_find_te test {NAME} : PASS\n"
	     "--------------------------------------------------------\n");
  
  return true;

}}"""

    cname = "tangle_slide_input_test{NAME}_combined.c".format(NAME=x)
    print "\tcombining " + doc_name + " and " + fname + "-> " + cname;
    
    bpfun = boilerplatefunction.format(NAME=x,DOCSTRING=docstring)
    print bpfun

    c_file = open(cname,'w')
    c_file.write("/* tangle interior test {name} assembled c code */\n\n".format(name=x))
    c_file.write(pdcode_output)
    c_file.write(bpfun)
    c_file.close()

    cnames.append(cname)
    test_names.append("pdint_check_tslide_data_ok_and_find_te_test{NAME}()".format(NAME=x))

print "assembling " + str(cnames) + " into test_tangle_slide_input_auto.c"
totalf = open('test_tangle_slide_input_auto.c','w')

prefixboilerplate = r"""
/*

   test_tangle_slide_input_auto.c : This code is auto-generated by 'assemble_tangleregentest.py' and shouldn't be edited by hand. It contains tests for the gatekeeper code in pd_tangle_slide.


*/

#ifdef HAVE_CONFIG_H
  #include"config.h"
#endif

#ifdef HAVE_STDIO_H
   #include<stdio.h>
#endif

#ifdef HAVE_ASSERT_H
   #include<assert.h>
#endif

#ifdef HAVE_STRING_H
   #include<string.h>
#endif

#ifdef HAVE_STDINT_H
   #include<stdint.h>
#endif

#ifdef HAVE_STDLIB_H
   #include<stdlib.h>
#endif

#ifdef HAVE_STDBOOL_H
   #include<stdbool.h>
#endif

#include<plcTopology.h>

int PD_VERBOSE=15;

/* We need to include a prototype for the function we're testing, because
   it's not exposed in the header files. */

bool pdint_check_tslide_data_ok_and_find_te(pd_code_t *pd,pd_tangle_t *t,
					    pd_idx_t n,
					    pd_idx_t *overstrand_edges, 
					    pd_idx_t *border_faces,
					    pd_idx_t **tangle_slide_edges,
                        bool *overstrand_goes_OVER);
"""

actualtest = r"""
int main() {{

  printf("test_tangle_slide_input (%s)\n",PACKAGE_STRING);
  printf("--------------------------------------------------------\n"
	 "Unit tests for pdint_check_tslide_data_ok_and_find_te(). \n"
	 "========================================================\n");

  if ({tests_string}) {{

    printf("=======================================================\n");
    printf("test_tangle_slide_input:  PASS.\n");
    exit(0);

  }} else {{

    printf("=====================================================\n");
    printf("test_tangle_slide_input:  FAIL.\n");
    exit(1);

  }}

  return 0;

}}"""

# We now need to assemble the {tests_string}

ts = "&&".join(test_names)

totalf.write(prefixboilerplate)

for x in cnames:
    tempf = open(x,'r')
    totalf.write(tempf.read())
    tempf.close()
    totalf.write("\n\n")

totalf.write(actualtest.format(tests_string=ts))
totalf.close

print "wrote test_tangle_slide_input_auto.c"


   
                 

    

    



