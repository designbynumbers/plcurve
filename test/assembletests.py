import glob
import subprocess
import re
pdstor_files = glob.glob('*.pdstor')
print ("Test Assembly Script\n"
      "==========================\n"
      "Found pdstor files...\n")
print "\n".join("\t"+str(x) for x in pdstor_files)
print "\nRunning carbonize...\n"
for x in pdstor_files:
    try:
        carbonize_results = subprocess.check_output(["carbonizepd",x])
        xroot = re.search('(.+)\.pdstor',x);
        print "\t" + x + "\t->\t" + xroot.group(0) + ".c"
    except subprocess.CalledProcessError:
        print "carbonize failed on " + x + "\n"
        print carbonize_results
        sys.exit(1)

print "\nExtracting test names and test type..."
names = [ ]
types = [ ]
for x in pdstor_files:
    m = re.match('(.+)_test(.+)_(.+).pdstor',x)
    names.append(m.group(2))
    types.append(m.group(1))
    
names = list(set(names))
names.sort()
types = list(set(types))
types.sort()
print "\tfound test names " + ", ".join(names)
print "\tfound test types " + ", ".join(types)

print "\nGenerating combined test files ... "



for x in names:
    
    fname = "./{type}_test{name}.c".format(type=types[0],name=x)
    before_name = "./{type}_test{name}_before.c".format(type=types[0],name=x)
    after_name = "./{type}_test{name}_after.c".format(type=types[0],name=x)
    doc_name = "./{type}_test{name}_doc.c".format(type=types[0],name=x)

    # now we need to parse the after pdstor file a little bit to see how many outputs are provided

    pdafter_name = "./{type}_test{name}_after.pdstor".format(type=types[0],name=x)
    pdafter_file = open(pdafter_name,'r')

    for line in pdafter_file:
        match = re.search('nelts (\d+)/(\d+) .+',line)
        if match:
            outputs = match.group(1)

    pdafter_file.close()

    #print 'types[0] is ' + types[0] + 'and types[0] == \'r2\' is ' + str(types[0] == 'r2') + ' and outputs == ' + outputs + ' and outputs == \'1\' is ' + str(outputs == '1')

    boilerplatefunction = ""

    if (types[0] == 'r1'):

        boilerplatefunction = """
bool {TYPE}test{NAME}() {{
  
  printf("---------------------------------------\\n"
	 "{TYPE} test {NAME}\n"
	 "---------------------------------------\\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_{TYPE}_test{NAME}_before_0();
  printf("done (passes pd_ok)\\n");

  printf("creating after configuration...");
  pd_code_t *after = pd_create_{TYPE}_test{NAME}_after_0();
  printf("done (passes pd_ok)\\n");
  
  printf("performing {TYPE} at crossing CROSS...");
  pd_code_t *newpd = DO_THING_HERE;
  if (!pd_ok(newpd)) {{ 
    printf("fail (output pd does not pass pd_ok)\\n");
    return false;
  }}
  printf("pass (output passes pd_ok)\\n");

  printf("testing for isomorphism with after configuration...");
  if (!pd_isomorphic(newpd,after)) {{ 
    pd_printf("fail (output pd %PD ",newpd);
    pd_printf("is not isomorphism to expected \\"after\\" configuration %PD\\n",after);
    return false;
  }}

  printf("pass (output and after isomorphic)\\n");
  
  printf("housecleaning...");
  pd_code_free(&before);
  pd_code_free(&after);
  pd_code_free(&newpd);
  printf("done\\n");

  printf("---------------------------------------\\n"
	 "{TYPE} test {NAME} : PASS\\n"
	 "---------------------------------------\\n");
  
  return true;

}}"""

    elif (types[0] == 'r2' and outputs == '2'):

        boilerplatefunction = """
    bool {TYPE}test{NAME}() {{
  
  printf("---------------------------------------\\n"
	 "{TYPE} test {NAME}\\n"
	 "---------------------------------------\\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_{TYPE}_test{NAME}_before_0();
  printf("done (passes pd_ok)\\n");

  printf("creating after-0 configuration...");
  pd_code_t *after0 = pd_create_{TYPE}_test{NAME}_after_0();
  printf("done (passes pd_ok)\\n");

  printf("creating after-1 configuration...");
  pd_code_t *after1 = pd_create_{TYPE}_test{NAME}_after_1();
  printf("done (passes pd_ok)\\n");

  pd_code_t *after[2] = {{after0,after1}};
    
  printf("performing r2 at crossings CROSS and CROSS...\\n");

  pd_idx_t  noutpd;
  pd_code_t **newpd;
  pd_idx_t  cr[2] = {{ CROSS,CROSS }};

  pd_R2_bigon_elimination(before,cr,&noutpd,&newpd);

  if (noutpd != 2) {{ 
    printf("fail (did not generate 2 output pd codes)\\n");
    return false;
  }}

  pd_idx_t ct;

  for(ct=0;ct<noutpd;ct++) {{ 

    if (!pd_ok(newpd[ct])) {{ 
      printf("fail (output pd %d does not pass pd_ok)\\n",ct);
      return false;
    }}
    printf("\\t output pd %d passes pd_ok...pass\\n",ct);
    printf("\\t testing for isomorphism with after configuration %d...",ct);
    if (!pd_isomorphic(newpd[ct],after[ct])) {{ 
      pd_printf("fail (output pd %PD ",newpd[ct]);
      pd_printf("is not isomorphism to expected \\"after\\" configuration %PD\\n",after[ct]);
      return false;
    }}

    printf("pass (output %d and after %d isomorphic)\\n",ct,ct);

  }}

  printf("performing r2 at crossings 4 and 5...pass\\n");
  printf("housecleaning...");

  pd_code_free(&before);
  pd_code_free(&after0);
  pd_code_free(&after1);
  pd_code_free(&newpd[0]);
  pd_code_free(&newpd[1]);

  free(newpd);

  printf("done\\n");

  printf("---------------------------------------\\n"
	 "{TYPE} test {NAME} : PASS\\n"
	 "---------------------------------------\\n");
  
  return true;

}}"""
        
    elif (types[0] == 'r2' and outputs == '1'):

        boilerplatefunction = """
    bool {TYPE}test{NAME}() {{
  
  printf("---------------------------------------\\n"
	 "{TYPE} test {NAME}\\n"
	 "---------------------------------------\\n");

  printf("creating before configuration...");
  pd_code_t *before = pd_create_{TYPE}_test{NAME}_before_0();
  printf("done (passes pd_ok)\\n");

  printf("creating after-0 configuration...");
  pd_code_t *after0 = pd_create_{TYPE}_test{NAME}_after_0();
  printf("done (passes pd_ok)\\n");

  pd_code_t *after[1] = {{after0}};
    
  printf("performing r2 at crossings CROSS and CROSS...\\n");

  pd_idx_t  noutpd;
  pd_code_t **newpd;
  pd_idx_t  cr[2] = {{CROSS,CROSS}};

  pd_R2_bigon_elimination(before,cr,&noutpd,&newpd);

  if (noutpd != 1) {{ 
    printf("fail (did not generate 1 output pd codes)\\n");
    return false;
  }}

  pd_idx_t ct;

  for(ct=0;ct<noutpd;ct++) {{ 

    if (!pd_ok(newpd[ct])) {{ 
      printf("fail (output pd %d does not pass pd_ok)\\n",ct);
      return false;
    }}
    printf("\\t output pd %d passes pd_ok...pass\\n",ct);
    printf("\\t testing for isomorphism with after configuration %d...",ct);
    if (!pd_isomorphic(newpd[ct],after[ct])) {{ 
      pd_printf("fail (output pd %PD ",newpd[ct]);
      pd_printf("is not isomorphism to expected \\"after\\" configuration %PD\\n",after[ct]);
      return false;
    }}

    printf("pass (output %d and after %d isomorphic)\\n",ct,ct);

  }}

  printf("performing r2 at crossings 4 and 5...pass\\n");
  printf("housecleaning...");

  pd_code_free(&before);
  pd_code_free(&after0);
  pd_code_free(&newpd[0]);

  free(newpd);

  printf("done\\n");

  printf("---------------------------------------\\n"
	 "{TYPE} test {NAME} : PASS\\n"
	 "---------------------------------------\\n");
  
  return true;

}}"""
    
    print "\tcombining " + before_name + ", " + after_name + ", (" + outputs +" outputs) and " + doc_name + " -> " + fname;
    bpfun = boilerplatefunction.format(TYPE=types[0],NAME=x)
    #print bpfun

    doc_file = open(doc_name,'r')
    before_file = open(before_name,'r')
    after_file = open(after_name,'r')
    f_file = open(fname,'w')

    f_file.write("/* {type} test {name} assembled c code */\n\n".format(type=types[0],name=x))

    f_file.write(before_file.read())
    f_file.write("\n\n")

    f_file.write(after_file.read())
    f_file.write("\n\n")

    f_file.write(doc_file.read())
    f_file.write("\n\n")

    f_file.write(bpfun)
    f_file.close()
                 

    

    



