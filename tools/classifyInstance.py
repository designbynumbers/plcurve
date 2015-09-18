#!/usr/bin/python
import argparse
import tempfile
import os
import subprocess
import re
import sys
import atexit
import time
import shutil

def garbagecollection(tdir):
    print "<in a real execution, would delete " + tdir + " now>"

suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
def humansize(nbytes):
    if nbytes == 0: return '0 B'
    i = 0
    while nbytes >= 1024 and i < len(suffixes)-1:
        nbytes /= 1024.
        i += 1
    f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
    return '%s %s' % (f, suffixes[i])

if __name__ == "__main__":
    

    parser = argparse.ArgumentParser(description='Run a piece of the larger classification job')
    parser.add_argument(
        'piecename',   
        help='.pdstor file to classify')    # by default this is a string argument; it should be
                                            # in the form "/inputs/<n>/<n>b-ddddd.pdstor", which
                                            # is the location of the file in the s3 bucket.

    args = parser.parse_args()

    # task 0 is to check that the input matches this form.

    if (re.match('/inputs/\d/\db-\d\d\d\d\d.pdstor',args.piecename) == None):

        print "cI: input filename " + (args.piecename) + " doesn't match format"
        sys.exit(1);

    print "cI: received valid input file " + args.piecename;

    # task 1 is to create a temp directory in which to work.

    print "cI: creating temp directory...",

    tempdir = tempfile.mkdtemp(prefix="classifierwork")
    os.chdir(tempdir)

    print "created " + tempdir;

    atexit.register(garbagecollection,tempdir)
    
    # task 2 is to download the file from the s3 bucket "knotclassification"

    print "cI: pulling input file from s3 bucket knotclassification...",

    try:
        awsOutput = subprocess.check_output("aws s3 cp s3://knotclassification" + args.piecename + " .",
                                            stderr=subprocess.STDOUT,
                                            shell=True)
    except subprocess.CalledProcessError as e:
        print "FAIL (s3 error, probably that the file doesn't exist)"
        print "\tinput: " + e.cmd
        print "\toutput: " + e.output
        sys.exit(1)

    inputpd = re.match('download: .+ to (.+)',awsOutput).group(1)
        
    try:
        size = os.stat("./" + inputpd).st_size
    except:
        print "FAIL (couldn't stat ./" + inputpd
        sys.exit(1)

    print "done (" + humansize(size) + ") input file ready"

    # task 3 is to actually run classify.py on the given input

    print "cI: running pdclassify on " + inputpd + "...",

    try:
        classifyOutput = subprocess.check_output("pdclassify " + inputpd,
                                                 stderr=subprocess.STDOUT,
                                                 shell=True)
    except subprocess.CalledProcessError as e:
        print "FAIL (pdclassify error)"
        print "\tinput: " + e.cmd
        print "\toutput: " + e.output
        sys.exit(1)

    ambigFile = inputpd + ".ambig.pdstor"
    countFile = inputpd + ".counts.tsv"
        
    try:
        ambigSize = os.stat(ambigFile).st_size
    except:
        print "FAIL (couldn't stat " + ambigFile
        sys.exit(1)

    try:
        countsSize = os.stat(countFile).st_size
    except:
        print "FAIL (couldn't stat " + countFile
        sys.exit(1)
        
    print "done (produced "+humansize(ambigSize)+" file of ambiguous knots, "+humansize(countsSize)+" file of counts)"

    # task 4 is to push the results into the s3 bucket

    print "cI: pushing results to s3 bucket knotclassification ...",

    basename = re.match("./(.+).pdstor",inputpd).group(1);
    n = re.match("./(\d+)b-.+.pdstor",inputpd).group(1);

    try:
        awsOutput = subprocess.check_output("aws s3 cp " + ambigFile + " s3://knotclassification/outputs/" + n + "/ambigs/" + basename + ".ambig.pdstor",
                                            stderr=subprocess.STDOUT,
                                            shell=True)
    except subprocess.CalledProcessError as e:
        print "FAIL (s3 error, probably couldn't find file to upload)"
        print "\tinput: " + e.cmd
        print "\toutput: " + e.output
        sys.exit(1)

    try:
        awsOutput = subprocess.check_output("aws s3 cp " + countFile + " s3://knotclassification/outputs/" + n + "/counts/" + basename + ".counts.tsv",
                                            stderr=subprocess.STDOUT,
                                            shell=True)
    except subprocess.CalledProcessError as e:
        print "FAIL (s3 error, probably couldn't find file to upload)"
        print "\tinput: " + e.cmd
        print "\toutput: " + e.output
        sys.exit(1)
        
    time.sleep(2) # wait 2 seconds for aws to realize that the upload is complete.

    # task 5 is to confirm that the upload actually happened
    
    try:
        awsOutput = subprocess.check_output("aws s3 ls s3://knotclassification/outputs/" + n + "/ambigs/"+ basename + ".ambig.pdstor",
                                            stderr=subprocess.STDOUT,
                                            shell=True)
        
    except subprocess.CalledProcessError as e:
        print "FAIL (s3 error, probably couldn't find file to uploadbucket)"
        print "\tinput: " + e.cmd
        print "\toutput: " + e.output
        sys.exit(1)

    match = re.match('.+ (\d+) '+ basename +'.ambig.pdstor',awsOutput)

    if (match == None) :
        print "FAIL (s3 couldn't confirm the file was uploaded)"
        print "couldn't match" + '.+ (\d)+ '+ basename +'.ambig.pdstor'
        print "awsOutput: " + awsOutput
        sys.exit(1)

    ambigCheckSize = int(match.group(1))

    if (int(ambigCheckSize) != int(ambigSize)):
        print "WARN (s3 uploaded size " + ambigCheckSize + " != local filesize " + ambigSize
        print "cI: This could be a size printing convention problem for large files"

    try:
        awsOutput = subprocess.check_output("aws s3 ls s3://knotclassification/outputs/" + n + "/counts/"+ basename + ".counts.tsv",
                                            stderr=subprocess.STDOUT,
                                            shell=True)
        
    except subprocess.CalledProcessError as e:
        print "FAIL (s3 error, probably couldn't find file to uploadbucket)"
        print "\tinput: " + e.cmd
        print "\toutput: " + e.output
        sys.exit(1)

    match = re.match('.+ (\d+) '+ basename +'.counts.tsv',awsOutput)

    if (match == None) :
        print "FAIL (s3 couldn't confirm the file was uploaded)"
        print "couldn't match" + '.+ (\d)+ '+ basename +'.ambig.pdstor'
        print "awsOutput: " + awsOutput
        sys.exit(1)

    countsCheckSize = int(match.group(1))

    if (int(countsCheckSize) != int(countsSize)):
        print "WARN (s3 uploaded size " + countsCheckSize + " != local filesize " + countsSize
        print "cI: This could be a size printing convention problem for large files"

        
    print "done (checked "+humansize(ambigCheckSize)+" file of ambigs uploaded,"+humansize(countsCheckSize)+" file of counts uploaded)"

    # task 6 is garbage collection (delete the temp directory)

    shutil.rmtree(tempdir)
    
