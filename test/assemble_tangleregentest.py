import glob
import subprocess
import re
import sys

pdstor_files = glob.glob('*.pdstor')
print ("Tangle Regeneration Test Assembly Script\n"
      "==========================\n"
      "Found pdstor files...\n")
print "\n".join("\t"+str(x) for x in pdstor_files)
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
    m = re.match('tangle_test(.+).pdstor',x)
    names.append(m.group(1))
    
names = list(set(names))
names.sort()

print "\tfound test names " + ", ".join(names)
print "\nGenerating combined individual test files ... "

cnames = []
test_names = []

for x in names:
    
    fname = "./tangle_test{name}.c".format(name=x)
    # the first thing we want to do is to slurp the contents of this file

    try:
        infile = open(fname,'r')
    except EnvironmentError:
        print "couldn't open " + fname + " for reading"
        sys.exit(1)
        
    pdcode_output = infile.read()
    infile.close()

    print "\tSlurped the c form of the pdstor from " + fname

    doc_name = "./tangle_test{name}.doc".format(name=x)

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
bool tangle_test{NAME}() {{

  {DOCSTRING}

  printf("---------------------------------------\n"
	 "tangle_regenerate test {NAME}\n"
	 "---------------------------------------\n");

  printf("creating pd...");
  pd_code_t *pd = pd_create_tangle_test{NAME}_0();
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

  printf("testing edge boundary orientations...");

  for(i=0;i<t->nedges;i++) {{

    if (t->edge_bdy_or[i] != edge_bdy_or[i]) {{

      pd_idx_t j;

      printf("fail (list of boundary orientations for tangles doesn't match expected at pos %d)\n",i);

      printf("found:   ");

      for(j=0;j<t->nedges;j++) {{

        if (t->edge_bdy_or[j] == in) {{ printf("  in "); }}
        else {{ printf(" out "); }}

      }}

      printf("\n");

      printf("expected:");

      for(j=0;j<t->nedges;j++) {{

        if (edge_bdy_or[j] == in) {{ printf("  in "); }}
        else {{ printf(" out "); }}

      }}

      printf("\n");

      return false;

    }}

  }}

  printf("pass (edge boundary orientations match expected)\n");

  printf("checking interior crossings...");
  qsort(interior_cross,(size_t)(ninterior_cross),sizeof(pd_idx_t),pd_idx_cmp);

  if (t->ninterior_cross != ninterior_cross) {{

    printf("fail (# interior crossings found (%d) != # expected (%d))\n",
    t->ninterior_cross,ninterior_cross);
    return false;

  }}

  for(i=0;i<t->ninterior_cross;i++) {{

     if (t->interior_cross[i] != interior_cross[i]) {{

        pd_idx_t j;

        printf("fail (list of interior crossings doesn't match expected at pos %d)\n",i);

        printf("found:   ");

        for(j=0;j<t->ninterior_cross;j++) {{

          printf(" %4d ",t->interior_cross[j]);

        }}

        printf("\n");

        printf("expected:");

        for(j=0;j<t->ninterior_cross;j++) {{

          printf(" %4d ",interior_cross[j]);

        }}

        printf("\n");

        return false;

    }}

  }}

  printf("pass (interior crossings match)\n");


  printf("checking interior edges...");
  qsort(interior_edge,(size_t)(ninterior_edges),sizeof(pd_idx_t),pd_idx_cmp);

  if (t->ninterior_edges != ninterior_edges) {{

    printf("fail (# interior edges found (%d) != # expected (%d))\n",
    t->ninterior_edges,ninterior_edges);
    return false;

  }}

  for(i=0;i<t->ninterior_edges;i++) {{

     if (t->interior_edge[i] != interior_edge[i]) {{

        pd_idx_t j;

        printf("fail (list of interior edges doesn't match expected at pos %d)\n",i);

        printf("found:    ");

        for(j=0;j<t->ninterior_edges;j++) {{

          printf(" %4d ",t->interior_edge[j]);

        }}

        printf("\n");

        printf("expected:");

        for(j=0;j<t->ninterior_edges;j++) {{

          printf(" %4d ",interior_edge[j]);

        }}

        printf("\n");

        return false;

    }}

  }}

  printf("pass (interior edges match)\n");
  
  printf("housecleaning...");
  
  pd_code_free(&pd);
  pd_tangle_free(&t);

  printf("done\n");

  printf("---------------------------------------\n"
         "tangle test {NAME} : PASS\n"
	     "---------------------------------------\n");
  
  return true;

}}"""

    cname = "tangle_test{NAME}_combined.c".format(NAME=x)
    print "\tcombining " + doc_name + " and " + fname + "-> " + cname;
    
    bpfun = boilerplatefunction.format(NAME=x,DOCSTRING=docstring)
    #print bpfun

    c_file = open(cname,'w')
    c_file.write("/* tangle interior test {name} assembled c code */\n\n".format(name=x))
    c_file.write(pdcode_output)
    c_file.write(bpfun)
    c_file.close()

    cnames.append(cname)
    test_names.append("tangle_test{NAME}()".format(NAME=x))

print "assembling " + str(cnames) + " into test_tangle_regenerate_auto.c"
totalf = open('test_tangle_regenerate_auto.c','w')

prefixboilerplate = r"""
/*

   test_tangle_regenerate_auto.c : This code is auto-generated by 'assemble_tangleregentest.py' and shouldn't be edited by hand. It contains tests for the orientation and interior code in pd_regenerate_tangle.


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

int PD_VERBOSE=50;
"""

actualtest = r"""
int main() {{

  printf("test_tangle_regenerate_auto (%s)\n",PACKAGE_STRING);
  printf("--------------------------------------------------------\n"
	 "Unit tests for pd_tangle_regenerate. \n"
	 "========================================================\n");

  if ({tests_string}) {{

    printf("=======================================================\n");
    printf("test_tangle_regenerate_auto:  PASS.\n");
    exit(1);

  }} else {{

    printf("=====================================================\n");
    printf("test_tangle_regenerate_auto:  FAIL.\n");
    exit(0);

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

print "wrote test_tangle_regenerate_auto.c"


   
                 

    

    



