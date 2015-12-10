#
# Master file for mkarxiv scripts.
# 

#
# Step 0: Paper-specific data.
#

my $stagingdir = "./arxiv/";
my $mainfile = "knotprobabilities_paper.tex";
my $newfile = "knot_probabilities_random_diagrams.tex";
#my @datafiles = qw(../TightKnot/best/ropelength_data.txt);  # Array of original filenames for data files (include full paths)

# Step 1 : Directory maintenance

use File::Path;
use File::Copy;

print("Clearing old directory...\n");
if (-d $stagingdir) {
  rmtree($stagingdir) or die("Couldn't clear old $stagingdir directory.");
}
mkdir($stagingdir) or die("Couldn't make new $stagingdir");

# Step 2 : Gather dependencies from TeX file and copy them to staging dir

print("Scanning $mainfile to gather names of required files....");
open(TEXFILE,$mainfile) or die("Couldn't open $mainfile\n");

my @inputs;
my @figs;
my $line;

while (<TEXFILE>) {

  $line = $_;
  if ($line =~ /\\input\{(.+?)\}/) { push(@inputs,$1); }
  if ($line =~ /\\includegraphics\[.+\]\{(.+?)\}/ or $line =~ /\\includegraphics\{(.+?)\}/) { push(@figs,$1); }
  if ($line =~ /.+\\includegraphics\[.+\]\{(.+?)\}/) { push(@figs,$1); }
  if ($line =~ /\\begin\{overpic\}\[.+\]\{(.+?)\}/) { push(@figs,$1); }

}

close(TEXFILE);
print("done.\n");

print ("File seems to depend on \\inputs \n\t");
print join("\n\t",@inputs);
print ("\n\nFile seems to depend on figures \n\t");
print join("\n\t",@figs);
print ("\n\nCopying these files to $stagingdir ... ");

foreach (@inputs) {

  copy($_,$stagingdir) or die("Couldn't copy $_ to $stagingdir\n");

}

#copy("diagrams.sty",$stagingdir) or die("Couldn't copy diagrams.sty to $stagingdir");

foreach (@figs) {

  if ($_ !~ /pdf/ && $_ !~ /png/ && $_ !~ /jpg/) { 

    # We go in order by extension:
    
    my $pdfname; 
    my $pngname;
    my $jpgname;
    
    $pdfname = "../../figs/$_.pdf";
    $pngname = "../../figs/$_.png";
    $jpgname = "../../figs/$_.jpg";
    
    copy($jpgname,$stagingdir) or copy($pdfname,$stagingdir) or 
      copy($pngname,$stagingdir) or die("Couldn't copy $pdfname or $pngname to $stagingdir\n");
    
  } else {

    copy("../../figs/$_",$stagingdir) or die("Couldn't copy ../../figs/$_ to $stagingdir\n");

  }

}

print ("done\n");

#
# Step 3. Copy $mainfile and the bbl files, make edits to mainfile as needed.
#

print ("Copying and editing $mainfile -> $newfile...");

open(NEWFILE,">${stagingdir}${newfile}") or die("Couldn't open ${stagingdir}${newfile}");
open(TEXFILE,$mainfile) or die("Couldn't open $mainfile\n");

my $line;

while (<TEXFILE>) {

  $line = $_;

  if ($line =~ /\\graphicspath|\\eads|\\pacs|\\submitto/ ) { next; } # don't copy the graphicspath line
  $line =~ s/figs\///g; # delete the figs/ prefix where it occurs.
  # $line =~ s/iopart/amsart/g; # change the style to article
  $line =~ s/\\ack\{\}/\\noindent Acknowledgements\./g;

  print NEWFILE $line;

}
close(NEWFILE);
close(TEXFILE);

print("done\n");

my $bbl;
my $newbbl;

$bbl = $mainfile;
$bbl =~ s/\.tex/\.bbl/g;

$newbbl = $newfile;
$newbbl =~ s/\.tex/\.bbl/g;

print("Copying $bbl -> ${stagingdir}${newbbl}...");
copy($bbl,"${stagingdir}${newbbl}") or die("Couldn't copy $bbl -> ${stagingdir}${newbbl}\n");
print("done\n");

#
# Stage 4: Error checking
#

print "Now checking result of pdflatex call on the new files...";

chdir($stagingdir) or die("Couldn't chdir to $stagingdir");

system("pdflatex -halt-on-error $newfile $2>/dev/null") == 0 or 
  print "FAIL. pdflatex error is:\n\n".`pdflatex -halt-on-error $newfile | tail`;

print "PASS.\n\n";

my $newroot;
$newroot = $newfile;
$newroot =~ s/\.tex//g;

unlink("$newroot.aux") or die("Couldn't clean $newroot.aux");
unlink("$newroot.pdf") or die("Couldn't clean $newroot.pdf");
unlink("$newroot.log") or die("Couldn't clean $newroot.log");


chdir("..") or die("Couldn't chdir to ..");

#
# Stage 6: Copy ancillary data to the /anc subdirectory.
#

my $ancdir = "./arxiv/anc/";

print("Clearing old directory...\n");
if (-d $ancdir) {
  rmtree($ancdir) or die("Couldn't clear old $ancdir directory.");
}
mkdir($ancdir) or die("Couldn't make new $ancdir");

copy("/Users/cantarel/randomdiagram/data/README-arxiv.txt","${ancdir}/README.txt") or die("Couldn't copy to ${ancdir}");
copy("/Users/cantarel/randomdiagram/data/supplementaldata/knot-table.csv","${ancdir}/knot-table.csv") or die("Couldn't copy knot table to ${ancdir}");
copy("/Users/cantarel/randomdiagram/data/supplementaldata/knot-frequency-03.csv","${ancdir}/knot-frequency-03.csv") or die("Couldn't copy knot-frequency-03.csv to ${ancdir}");
copy("/Users/cantarel/randomdiagram/data/supplementaldata/knot-frequency-04.csv","${ancdir}/knot-frequency-04.csv") or die("Couldn't copy knot-frequency-04.csv to ${ancdir}");
copy("/Users/cantarel/randomdiagram/data/supplementaldata/knot-frequency-05.csv","${ancdir}/knot-frequency-05.csv") or die("Couldn't copy knot-frequency-05.csv to ${ancdir}");
copy("/Users/cantarel/randomdiagram/data/supplementaldata/knot-frequency-06.csv","${ancdir}/knot-frequency-06.csv") or die("Couldn't copy knot-frequency-06.csv to ${ancdir}");
copy("/Users/cantarel/randomdiagram/data/supplementaldata/knot-frequency-07.csv","${ancdir}/knot-frequency-07.csv") or die("Couldn't copy knot-frequency-07.csv to ${ancdir}");
copy("/Users/cantarel/randomdiagram/data/supplementaldata/knot-frequency-08.csv","${ancdir}/knot-frequency-08.csv") or die("Couldn't copy knot-frequency-08.csv to ${ancdir}");
copy("/Users/cantarel/randomdiagram/data/supplementaldata/knot-frequency-09.csv","${ancdir}/knot-frequency-09.csv") or die("Couldn't copy knot-frequency-09.csv to ${ancdir}");
copy("/Users/cantarel/randomdiagram/data/supplementaldata/knot-frequency-10.csv","${ancdir}/knot-frequency-10.csv") or die("Couldn't copy knot-frequency-10.csv to ${ancdir}");

#
# Stage 6: Package for arxiv upload
#

my $tarname;

$tarname = $newfile;
$tarname =~ s/\.tex/-arxiv.tgz/g;

print "Making $tarname for arxiv upload...";

system("tar -czf $tarname $stagingdir") == 0 or die("tar call failed");

print "done.\n";

my $tarsize;
$tarsize = `du -h $tarname`;
$tarsize =~ s/\s.+$//g;
chomp($tarsize);

print "Upload tarfile $tarname of size $tarsize created and ready for arxiv!\n";







  
       


