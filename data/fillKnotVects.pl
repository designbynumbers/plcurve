#!perl
#
#
# Fill the knotvects directory with files copied from TightKnot/vects in
# order to compute a table of HOMFLY polynomials for plc_knottype.
#

use strict;
use warnings;
use Lingua::EN::Numbers qw(num2en num2en_ordinal);
use File::Slurp;
use Data::Dumper;
use File::Copy;

sub make_safekey {

    my $key = $_[0];
    my $safekey;
    my $texkey;

    my $cr;
    my $ind;
    my $comp;

    my $safeone;
    my $safetwo;
    my $safethree;
    
    if ( $key !~ /.vect/ ) {  #this is a knot type
	
	if ( $key =~ /(\d+)_(\d+)_(\d+)/ ) { # a link
	    
	    $safeone = $1; $safetwo = $2; $safethree = $3;

	    $cr = num2en($safeone); $comp = num2en($safetwo); $ind = num2en($safethree);
	    $cr =~ s/\s/\-/g; $comp =~ s/\s/\-/g; $ind =~ s/\s/\-/g;
	   
	    $cr =~ s/\b(\w)(\w*)/uc($1).lc($2)/ge;
	    $comp =~ s/\b(\w)(\w*)/uc($1).lc($2)/ge;
	    $ind =~ s/\b(\w)(\w*)/uc($1).lc($2)/ge;
	    
	    $safekey = $cr.$comp.$ind;
	    $safekey =~ s/-//g;
	    $texkey = $safeone."^{".$safetwo."}_{".$safethree."}";

	    $cr = $safeone; $comp = $safetwo; $ind = $safethree;
	    
	} elsif ( $key =~ /(\d+)_(\d+)/ ) { # a knot
	
	    $safeone = $1; 
	    $safetwo = $2;
	    
	    $texkey = $safeone."_{".$safetwo."}";
    
	    $cr = num2en($safeone); $ind = num2en($safetwo);
	    $cr =~ s/\s/\-/g; $ind =~ s/\s/\-/g;

	    $cr =~ s/\b(\w)(\w*)/uc($1).lc($2)/ge;
	    $ind =~ s/\b(\w)(\w*)/uc($1).lc($2)/ge;
	    

	    $safekey = $cr.$ind;
	    $safekey =~ s/-//g;
	
	    $cr = $safeone; $comp = 1; $ind = $safetwo;

	} else {

	    print "FAILED to match key ".$key."\n";

	}


    }
	
    my @retstuff = ($safekey,$texkey,$cr,$comp,$ind);

    return @retstuff;
	    
}

#open(DB,'KnotDB.pl');
my $DBFILE = read_file('../../TightKnot/KnotDB.pl');
my %knots;

eval($DBFILE);

#close(DB);

my $numkeys = keys ( %knots );
print "We have $numkeys entries in the original database. Looking for knots.\n";

# first, a pass to remove files from the hash

my $key;

foreach $key ( sort(keys %knots) ) {

    if ($key =~ /.vect/) {

	delete $knots{$key};

    }

}

$numkeys = keys ( %knots );
print "There seem to be $numkeys knots and links in the database...\n";


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

foreach $knotline ( @rawtable ) {

    chomp($knotline);
    ($knot,$symtype) = split(/ /,$knotline);

    print("$knot -> ");

    # We now need to split the "m", "r", or "rm" from the knot type

    if ($knot =~ m/(\d+_\d+)(m|rm|r)/) {
	
	$tag = $1; $mod = $2;
	if ($mod =~ m/rm/) {

	  $operation = '"(-1,-1,e)"';

	} elsif ($mod =~ m/m/) {

	  $operation = '"(-1,1,e)"';

	} elsif ($mod =~ m/r/) {

	  $operation = '"(1,-1,e)"';

	}

    } else {

	$tag = $knot; $mod = 'e';
	$operation = '"(1,1,e)"';

    }

    print("$tag and $mod -> $operation");

    # Now we need to create the appropriate file and move it into the vect
    # directory.

    $fname = $knots{$tag}{filename};

    copy("/Users/cantarel/TightKnot/vects/".$fname,"./knotvects/".$tag.".vect") or die("Can't copy file $fname from TightKnot/vects/ to ./knotvects/$tag.vect");

    chdir("./knotvects");

    open(WHITTENVECT,"whittenvect ".'2>/dev/null'." -o $operation $tag.vect |");
    while (<WHITTENVECT>) {
      if (/Wrote modified VECT file to:\s*(\S+)\n/) { 
	    $newfile = $1;
	  } else {
	    $newfile = "Could not detect filename";
	  } 
    }
    close WHITTENVECT;

    chdir("..");

    print(" -> $newfile\n");

    # Now we need to actually generate the header file with the data by running the homfly code. 
	   
  }




