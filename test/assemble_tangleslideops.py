import glob
import subprocess
import re
import sys
import os.path

pdstor_files = glob.glob('*_before.pdstor')
print ("Tangle Slide Operation Test Assembly Script\n"
      "==========================\n"
      "Found pdstor files...\n")
print "\n".join("\t"+str(x) for x in pdstor_files)
print "Looking for corresponding .c files and moving to .doc..."
for x in pdstor_files:
    m = re.match('(.+)_before.pdstor',x)
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

pdstor_files = glob.glob('*.pdstor')                
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

pdstor_files = glob.glob('*_before.pdstor')                        
print "\nExtracting test names..."
names = [ ]
for x in pdstor_files:
    m = re.match('tangle_slide_operation_test(.+)_before.pdstor',x)
    names.append(m.group(1))
    
names = list(set(names))
names.sort()

print "\tfound test names " + ", ".join(names)
print "\nGenerating combined individual test files ... "

cnames = []
test_names = []

for x in names:
    
    fname_before = "./tangle_slide_operation_test{name}_before.c".format(name=x)
    # slurp the contents of this file

    try:
        infile = open(fname_before,'r')
    except EnvironmentError:
        print "couldn't open " + fname_before + " for reading"
        sys.exit(1)
        
    pdcode_before = infile.read()
    infile.close()

    print "\tSlurped c form of the pdstor from " + fname_before

    fname_after = "./tangle_slide_operation_test{name}_after.c".format(name=x)
    # slurp the contents of this file

    try:
        infile = open(fname_after,'r')
    except EnvironmentError:
        print "couldn't open " + fname_after + " for reading"
        sys.exit(1)
        
    pdcode_after = infile.read()
    infile.close()

    print "\tSlurped c form of the output pdstor from " + fname_after

    children = []
    print "parsing output pdstor to find nchildren..."
    for m in re.finditer("pd_code_t \*(.+)\(\) {",pdcode_after):
        children.append(m.group(1));
    nchildren = len(children)
    print "\tFound " + str(nchildren) + " children " + ",".join(children)
    buildchildrenstring = ["pd_idx_t nxchildren = " + str(nchildren) + ";\n","  pd_code_t *xchildren[" + str(nchildren) + "];\n"];

    count =0;
    for c in children:
        buildchildrenstring.append("  xchildren["+str(count)+"] = " + c + "();\n");
        count += 1;

    xchildren = " ".join(buildchildrenstring)
        
    doc_name = "./tangle_slide_operation_test{name}.doc".format(name=x)

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
bool pdint_check_tangle_slide_ops_test{NAME}() {{

  {DOCSTRING}

  printf("--------------------------------------------------\n"
    	 "pdint check tangle slide ops test {NAME}\n"
         "--------------------------------------------------\n");

  printf("creating 'before' pd...");
  pd_code_t *pd = pd_create_tangle_slide_operation_test{NAME}_before_0();
  if (!pd_ok(pd)) {{

    printf("fail (doesn't pass pd_ok)\n");
    return false;

  }}
  
  printf("pass (passes pd_ok)\n");

  printf("parsing tangle information...");
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
  printf("parsing overstrand_edges and border_faces...");
  printf("%d overstrand edges\n",noverstrand_edges);

  printf("running pd_tangle_slide...");

  pd_code_t **children;
  pd_idx_t nchildren;

  pd_tangle_slide(pd,t,
                  noverstrand_edges,overstrand_edges,border_faces,
                  &nchildren,&children);

  printf("pass (didn't crash, returned %d children)\n",nchildren);
  
  printf("checking children for pd_ok...\n");
  
  for(i=0;i<nchildren;i++) {{

     if (!pd_ok(children[i])) {{

        pd_printf("child %d fails pd_ok. \n %PD",children[i],i);
        return false;

     }} else {{

        printf("\t%d crossing child pd #%d is pd_ok...pass\n",
               children[i]->ncross,i);

     }}

   }}
  
   {LOADXCHILDREN}

   printf("checking %d expected children for pd_ok...\n",nxchildren);

   for(i=0;i<nxchildren;i++) {{

    if (!pd_ok(xchildren[i])) {{

       pd_printf("expected child %d failed pd_ok:\n %PD",xchildren[i],i);
       return false;

    }} else {{

       printf("\t%d crossing expected child pd #%d is pd_ok...pass\n",
              xchildren[i]->ncross,i);

    }}

   }}

   printf("comparing list of xchildren against actual children...\n");

   if (!compare_list_of_pds(nchildren,children,nxchildren,xchildren)) {{

       printf("FAIL. Couldn't match expected and actual children\n"
              "after tangle_slide\n");
       return false;

   }} else {{

       printf("comparing list of xchildren against actual children...pass\n");

    }}

  printf("housecleaning...\n");

  printf("\tfreeing parent pd and tangle...");
  pd_code_free(&pd);
  pd_tangle_free(&t);
  printf("done\n");

  printf("\tfreeing list of children and children buffer...");
  for(i=0;i<nchildren;i++) {{
     pd_code_free(&(children[i]));
  }}
  free(children);
  printf("done\n");

  printf("\tfreeing list of xchildren...");
  for(i=0;i<nxchildren;i++) {{
     pd_code_free(&(xchildren[i]));
  }}
  printf("done\n");

  printf("--------------------------------------------------------\n"
         "pdint check tangle slide ops test{NAME} : PASS\n"
	     "--------------------------------------------------------\n");
  
  return true;

}}"""

    cname = "tangle_slide_ops_test{NAME}_combined.c".format(NAME=x)
    print "\tcombining " + doc_name + ", " + fname_before + " and " + fname_after + "-> " + cname;
    
    bpfun = boilerplatefunction.format(NAME=x,DOCSTRING=docstring,LOADXCHILDREN=xchildren)
    print bpfun

    c_file = open(cname,'w')
    c_file.write("/* tangle slide operations {name} assembled c code */\n\n".format(name=x))
    c_file.write(pdcode_before)
    c_file.write(pdcode_after)
    c_file.write(bpfun)
    c_file.close()

    cnames.append(cname)
    test_names.append("pdint_check_tangle_slide_ops_test{NAME}()".format(NAME=x))

print "assembling " + str(cnames) + " into test_tangle_slide_ops_auto.c"
totalf = open('test_tangle_slide_ops_auto.c','w')

prefixboilerplate = r"""
/*

   test_tangle_slide_ops_auto.c : This code is auto-generated
   by 'assemble_tangleslideops.py' and shouldn't be edited
   by hand. It contains tests for pd_tangle_slide.


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
#include<pd_multidx.h>
#include<pd_perm.h>

extern int PD_VERBOSE;

/* We need a utility function to do the tests... */

bool compare_list_of_pds(pd_idx_t nA, pd_code_t **A,
                         pd_idx_t nB, pd_code_t **B)

/* Figure out whether there's a match between the list of
   pds in A and the list in B. Since we expect the lists
   to be kind of small, we're just brute-forcing
   this (of course, there are many cleverer ways to do it). */

{{
  printf("\tchecking that lists of pd codes are same size...");
  
  if (nA == nB) {{ 

    printf("pass (%d child pd codes)\n",nA);

  }} else {{ 

    printf("FAIL. (lists are size %d != %d)\n",nA,nB);
    return false;

  }}

  pd_perm_t *perm;
  perm = pd_new_perm((void *)(&nA));
  bool all_iso = true;

  printf("\tchecking up to %u perms for diagram isotopy match between lists...",
	 pd_nperms((void *)(perm)));
  
  do {{  /* Check the current permutation. */

    pd_idx_t i; 

    for(i=0,all_iso = true;i < nA && all_iso;i++) {{

      if (!pd_diagram_isotopic(A[i],B[perm->map[i]])) {{ all_iso = false; }}

    }}
    
    pd_increment_perm((void *)(perm));

  }} while (!all_iso && !pd_perm_is_e(perm));

  if (all_iso) {{ 

    /* We need a trick here-- we can't DECREMENT the perm
       (though we'd like to) in order to report the
       permutation that we actually used. 
       But we can increment it n-1 times... */

    pd_idx_t np = pd_nperms((void *)(perm));
    pd_idx_t i;
    for(i=0;i<np-1;i++) {{ pd_increment_perm((void *)(perm)); }}

    pd_printf("pass\n\t(%PERM matches pd codes "
              "in list A with list B)\n",NULL,perm);
    
  }} else {{

    printf("FAIL (could not find a match between codes in list)\n");
    return false;

  }}

  pd_free_perm((void **)(&perm));

  return true;

}}                        
"""

actualtest = r"""
int main() {{

  PD_VERBOSE = 15;

  printf("test_tangle_slide_operation (%s)\n",PACKAGE_STRING);
  printf("--------------------------------------------------------\n"
	 "Unit tests for pd_tangle_slide. \n"
	 "========================================================\n");

  if ({tests_string}) {{

    printf("=======================================================\n");
    printf("test_tangle_slide_operation:  PASS.\n");
    exit(0);

  }} else {{

    printf("=====================================================\n");
    printf("test_tangle_slide_operation:  FAIL.\n");
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

print "wrote test_tangle_slide_ops_auto.c"


   
                 

    

    



