#!perl
#
#
# make_homfly_table.pl : Creates a table of homfly polynomials for prime and composite links for plCurve to use in classifying curves by knot type.
#
#

use warnings;
use Lingua::EN::Numbers qw(num2en num2en_ordinal);
use File::Slurp;
use Data::Dumper;
use File::Copy;


# now a pass through the knotchart-prime.txt file.
my @rawtable;

open(TABLE,"knotchart-prime.txt") || die("Couldn't open knotchart...\n");
@rawtable = <TABLE>;
close(TABLE);

my $knotline;
my $knot;
my $tag;
my $symtype;
my $mod;
my $operation;
my $newfile;
my $fname;
my @ktab;
my $homfly;

foreach $knotline ( @rawtable ) {

    my $knotrec = { };
    my $storedop;

    chomp($knotline);
    ($knot,$symtype) = split(/ /,$knotline);

    print("$knot -> ");

    # We now need to split the "m", "r", or "rm" from the knot type

    if ($knot =~ m/(\d+_\d+)(m|rm|r)/) {
	
	$tag = $1; $mod = $2;
	if ($mod =~ m/rm/) {

	  $operation = "[-1,-1,e]";
	  $storedop = '"(-1,-1,e)"';

	} elsif ($mod =~ m/m/) {

	  $operation = "[-1,1,e]";
	  $storedop = '"(-1,1,e)"';

	} elsif ($mod =~ m/r/) {

	  $operation = "[1,-1,e]";
	  $storedop = '"(1,-1,e)"';

	}

    } else {

	$tag = $knot; $mod = 'e';
	$operation = "[1,1,e]";
	$storedop = '"(1,1,e)"';

    }

    print("$tag and $mod -> $operation");

    # Now we need to identify the corresponding file in ./knotvects/

    $fname = "$tag.$operation.vect";

    print (" -> $fname");

    chdir("./knotvects");

    unless (-f $fname) { die("Couldn't locate prime file $fname in ./knotvects/");}

    # Now we need to run the homfly code. 

    open(KNOTTYPE,"knottype -qh /Users/cantarel/plCurve/data/knotvects/$fname |");
    while (<KNOTTYPE>) {

      if (/Homfly polynomial:\((.+)\)/) {
	$homfly = $1;
      } 
    
    }

    print(" -> $homfly\n");
	   
    $knotrec->{homfly} = $homfly;
    @{$knotrec->{sym}}    = ($storedop);

    $tag =~ /(\d+)_(\d+)/;
    @{$knotrec->{cr}} = ($1);
    @{$knotrec->{ind}} = ($2);
    @{$knotrec->{filename}} = ($fname);
    $knotrec->{nsummands} = 1;

    push(@ktab,$knotrec);

   # print Dumper(%{$knotrec})."\n";

    chdir("..");

  }

# We now generate a table of composites. 

open(TABLE,"knotchart-composite-only.txt") || die("Couldn't open knotchart...\n");
@rawtable = <TABLE>;
close(TABLE);

foreach $knotline (@rawtable) {

   my $knotrec = { };
   my @knotline;

   chomp($knotline);
   @knotline=split(/\s/,$knotline);
   $symtype = pop(@knotline);
   # clean up the knotline
   my @hold = ( );
   foreach ( @knotline ) { 
     if (m/(\d+)_(\d+)/) {push(@hold,$_);} 
   }
   @knotline = @hold;
 
   my $summand;
   my $storedop;
   my @files = ( );
   my @cr = ( );
   my @ind = ( );
   my @sym = ( );

   print scalar @knotline." summands ->";

   foreach $summand (@knotline) {

     if ($summand =~ m/(\d+_\d+)(m|rm|r)/) {
	
	$tag = $1; $mod = $2;
	if ($mod =~ m/rm/) {

	  $operation = "[-1,-1,e]";
	  $storedop = '"(-1,-1,e)"';

	} elsif ($mod =~ m/m/) {

	  $operation = "[-1,1,e]";
	  $storedop = '"(-1,1,e)"';

	} elsif ($mod =~ m/r/) {

	  $operation = "[1,-1,e]";
	  $storedop = '"(1,-1,e)"';

	}

      } else {
	
	$tag = $summand; $mod = 'e';
	$operation = "[1,1,e]";
	$storedop = '"(1,1,e)"';

      }
   
     push(@files,'"'."/Users/cantarel/plCurve/data/knotvects/$tag.$operation.vect".'"');
     
     $tag =~ m/(\d+)_(\d+)/;

     push(@cr,$1);
     push(@ind,$2);
     push(@sym,$storedop);

   }

   print "@knotline -> homfly ->";

   # We now compute the HOMFLY for the composite

   $homfly = "Error";

   open(KNOTTYPE,"knottype -q#h @files |");
   while (<KNOTTYPE>) {
     
      if (/Composite Homfly polynomial:\((.+)\)/) {
	$homfly = $1;
      } 

      if (/open/) {
	print;
	die;
      }
      
    }

   close (KNOTTYPE);

   if ($homfly eq "Error" && $cr[0] eq "3" && $cr[1] eq "4" && $ind[0] eq "1" && $ind[1] eq "1") {

     $homfly = "[[2] 0 3 0 3 0 1]N[0]N[-1] 0 -3 0 -2N[0]N[0] 0 1N ";

   }

    print(" $homfly\n");

   # Now we push everything into the hash and into the knot table.

     $knotrec->{homfly}    = $homfly;
   @{$knotrec->{sym}}      = @sym; 
   @{$knotrec->{cr}}       = @cr;
   @{$knotrec->{ind}}      = @ind;
   @{$knotrec->{filename}} = @files;
     $knotrec->{nsummands} = scalar @knotline;
   
   push(@ktab,$knotrec);

  # print Dumper(%{$knotrec})."\n";
   
   chdir("..");
   
   }

print "Completed processing knotchart-prime and knotchart-composite-only...\n";
print "Adding unknot...\n";

$knotrec->{homfly} = "[[0]]N";
@{$knotrec->{sym}} = ('"(1,1,e)"');
@{$knotrec->{cr}}  = ("0");
@{$knotrec->{ind}} = ("1");
@{$knotrec->{filename}} = ("None");
$knotrec->{nsummands} = 1;

push(@ktab,$knotrec);

# We have now generated a complete list of homflies for everything in the prime and composite tables. 
# Our next mission (should we choose to accept it!) is to output this stuff to a header file. 

print "Completed adding unknot...\n";

# Now we sort the table by homfly

my @knotlist = sort{ $a->{homfly} cmp $b->{homfly} } @ktab ;
my $ntypes;

$ntypes = scalar @ktab;
$ntypes++;

print "Generated table of $ntypes homfly polynomials for different knot types. Writing header file...\n";

# We generate the homfly database file

open(HOMFLY,">/Users/cantarel/plCurve/homfly.h") || die("Could not open homfly.h header file.\n");

print HOMFLY <<ENDC;

/* 
    homfly.h : This file contains a small data structure designed to store the homfly
    polynomials of knots and links to 10 crossings. This file is automatically generated
    by /data/make_homfly_table.pl and should not be edited by hand. We are filling in values
    of data type plc_knottype, which is defined in plCurve.h.

    #define MAXPRIMEFACTORS 10
    #define MAXHOMFLY       1024

    struct knottypestruct {

      int  nf;                            // Number of prime factors 
      int  cr[MAXPRIMEFACTORS];           // Crossing number of each prime factor 
      int  ind[MAXPRIMEFACTORS];          // Index (in Rolfsen or Cerf) of each prime factor 
      char sym[MAXPRIMEFACTORS][128];     // Symmetry tag (Whitten group element) for each prime factor 
      char homfly[MAXHOMFLY];             // Homfly polynomial (as plc_lmpoly output) 

    } plc_knottype;

*/

ENDC

print HOMFLY "#define KTDBSIZE ".($ntypes-1)."\n";

print HOMFLY "plc_knottype ktdb[$ntypes]  = { \\\n";

my %tag;

foreach $tag ( @knotlist ) {

 # print "@{$tag->{cr}}\n";
 # print Dumper($tag)."\n";
 # print join(',',@{$tag->{cr}})."\n";

  print HOMFLY "{ $tag->{nsummands}, { ".join(',',@{$tag->{cr}})."}, { ".join(',',@{$tag->{ind}})." }, { ".join(',',@{$tag->{sym}})." }, \"$tag->{homfly}\" }, \\\n";

}

print HOMFLY " { 1, { 0 }, { 1 }, { \"(1,1,e)\"}, {\"[[0]]N\"} } \\\n";  # add a guy to take care of trailing comma
print HOMFLY "};\n";

close(HOMFLY);

   

     


