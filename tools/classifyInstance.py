import argparse
import tempfile
import os
import subprocess
import re
import sys
import atexit


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

    print "done"
    print "cI: pdclassify output\n" + classifyOutput


    
    
    
    
