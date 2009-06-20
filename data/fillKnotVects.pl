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

open(TABLE,"knotchart-prime.txt") || die("Couldn't open knotchart...\n");
@rawtable = <TABLE>;
close(TABLE);

my $knotline;
my $knot;
my $tag;
my $symtype;
my $mod;

foreach $knotline ( @rawtable ) {

    chomp($knotline);
    ($knot,$symtype) = split(/ /,$knotline);

    # We now need to split the "m", "r", or "rm" from the knot type

    if ($knot =~ /(\d+_\d+)(\c+)/) {
	
	$tag = $1; $mod = $2;

    } else {

	$tag = $knot; $mod = 'e';

    }

    # Now we need to create the appropriate file and move it into the vect
    # directory.

    


# now a pass to generate files 

my $safekey;
my $texkey;
my $cr;
my $comp;
my $ind;

foreach $key ( sort(keys %knots) ) {

    ($safekey,$texkey,$cr,$comp,$ind) = &make_safekey($key);

    # Now we copy and rename files as needed. 

    if ($comp == 1) {

	# We now load the symmetry data from knotchart-prime.txt

	

    $knots{$key}{ropelength} = sprintf("%2.4f",$knots{$key}{ropelength} + 0.00005);
    $knots{$key}{upperbound} = sprintf("%2.4f",$knots{$key}{upperbound} + 0.00005);

    print TEX "\\newcommand{\\".$safekey."}{".$texkey."}\n";
    print TEX "\\newcommand{\\".$safekey."RRUB}{".$knots{$key}{upperbound}."}\n";
    print TEX "\\newcommand{\\".$safekey."RProp}{".$knots{$key}{ropelength}."}\n";
    
    $knots{$key}{tex} = $texkey;
    $knots{$key}{cr} = $cr;
    $knots{$key}{comp} = $comp;
    $knots{$key}{index} = $ind;

    print "Assigned key to ".$key." of ".$texkey.", cr of ".$cr.", comp of ".$comp.", index of ".$ind." with Latex key of ".$safekey."\n";
   
}

close TEX;

# We first delete manually some troublesome cases

delete $knots{'0_1'};
delete $knots{'0_2_1'};
delete $knots{'0_3_1'};
delete $knots{'0_4_1'};
delete $knots{'3_1_1'};

my @knotlist = sort{ 

		     if (($knots{$a}{cr} <=> $knots{$b}{cr}) eq 0 ) {

			 if (($knots{$a}{comp} <=> $knots{$b}{comp}) eq 0 ) {

			     return $knots{$a}{index} <=> $knots{$b}{index};

			 }

			 return $knots{$a}{comp} <=> $knots{$b}{comp};

		     }

		     return $knots{$a}{cr} <=> $knots{$b}{cr};

		     } keys( %knots ) ;

print "Now organizing the table by crossing number...\n";

my @bycrossing;

foreach $key ( @knotlist ) {

    if (defined $knots{$key}{cr}) { push(@{ $bycrossing[$knots{$key}{cr}] },$key); }

}

print "Example: 5 crossing knots are @{ $bycrossing[5] }. \n";
print "Example: 9 crossing knots are @{ $bycrossing[9] }. \n";

print "Now making list of knots and links....\n";

my @byknottype;
my $crnum;
my $lastkey;

foreach $crnum (2..10) {

    foreach $key ( @{ $bycrossing[$crnum] } ) {

	# If the # of components changes, output a trimmed midrule to highlight it.

	if (defined($lastkey)) {

	    if ($knots{$lastkey}{comp} != $knots{$key}{comp}) {

		#my $entry = pop(@byknottype);
		#chomp($entry);
		#push(@byknottype,$entry.' \midrule'."\n");

		push(@byknottype,'\addlinespace[0.45em]\cmidrule(r{1.5em}l{1.5em}){1-3}\addlinespace[0.45em]'."\n");

	    }

	}

	$lastkey = $key;

	# Now output the regular line. 

	push(@byknottype,'$'.$knots{$key}{tex}.'$ & $'.$knots{$key}{ropelength}.'$ & $'.$knots{$key}{upperbound}.'$ \\\\'."\n");

    }

    # When crossing number changes, insert an untrimmed midrule.

    $lastkey = undef;
#    my $entry = pop(@byknottype);
#    chomp($entry);
#    push(@byknottype,$entry.' \midrule'."\n");

    if ($crnum != 10) {push(@byknottype,'\addlinespace[0.45em]\cmidrule(r{0.5em}l{0.5em}){1-3}\addlinespace[0.45em]'."\n");}
    
}

my $colheads = 'Link & $\PRop$ & $\Rop$';
my @rows = qw(45 45 45 43);
my $caption = "Ropelengths of Tight Knots and Links by Knot Type";

# now a pass to generate the first actual table

&perltable("KnotTable.tex",\@rows,\@byknottype,$colheads,$caption,3,"ByKnotTable");    

# ############  Now we output the ordered by ropelength table.

@knotlist = sort{ $knots{$a}{upperbound} <=> $knots{$b}{upperbound} } keys( %knots );

my @byropelength;

foreach $key ( @knotlist ) {

    push(@byropelength,'$'.$knots{$key}{tex}.'$ \\\\'."\n");

}

@rows = qw(45 45 45);
$colheads = 'Link';
$caption = "Knot and Link Types sorted by Ropelength";

&perltable("KnotTableByRop.tex",\@rows,\@byropelength,$colheads,$caption,8,"ByRopTable");

#
# Now we output the strutplots
#

open(BESTIARY,">Bestiary.tex") or die("Couldn't open Bestiary.tex\n"); 
my $keynum = 0;

foreach $crnum (2..10) {
    
    print "Now outputting pagerefs for $crnum crossing knots and links...\n";
    
    print BESTIARY '\chapter{$'.$crnum.'$ crossing knots and links}'."\n";

    my $prcap = 'Page Numbers of $'.$crnum.'$ crossing knots and links';
    my @pagerefs = ();

    my @indexkeys;

    foreach my $key ( @{ $bycrossing[$crnum] } ) {

	($safekey,$texkey,$cr,$comp,$ind) = &make_safekey($key);		
	push(@pagerefs,'$'.$texkey.'$, p. \pageref{'.$safekey.'Page} \\\\'."\n");

    }
    
    my @prrows = qw(40 40 40 40 40);
    my $prcolhead = "Link";

    my $tfname = "Cr".$crnum."Table.tex";

    &perltable($tfname,\@prrows,\@pagerefs,$prcolhead,$prcap,5,"PageRefTab".$crnum);

    print BESTIARY '\input{'.$tfname.'}'."\n\n";

     # We have now output the chapter table. We go ahead and make pages now.

    print "    Writing knot ";
     
    foreach my $localkey ( @{ $bycrossing[$crnum] } ) {

	$keynum++;
	print " $localkey ";

	($safekey,$texkey,$cr,$comp,$ind) = &make_safekey($localkey);
	
	print BESTIARY '\newpage'."\n";
	print BESTIARY '\label{'.$safekey.'Page}'."\n";
	if ($comp == 1) { 

	    print BESTIARY '\index{'.$cr.'@'.$cr.' crossing links!'.$comp.'@with '.$comp.' component!'.$keynum.'@$'.$texkey.'$}'."\n";

	} else {

	    print BESTIARY '\index{'.$cr.'@'.$cr.' crossing links!'.$comp.'@with '.$comp.' components!'.$keynum.'@$'.$texkey.'$}'."\n";

	}

	# Now include the mugshots
        
	if (!defined $knots{$localkey}{mugshot}) { 

	    die("Mugshot not defined for $localkey. Data is".Dumper($knots{$localkey})); 

	}

	my $textvar = '}.jpg';

	$knots{$localkey}{mugshot}[0] =~ s/.jpg/$textvar/;
	$knots{$localkey}{mugshot}[1] =~ s/.jpg/$textvar/;
	$knots{$localkey}{mugshot}[2] =~ s/.jpg/$textvar/;

	print BESTIARY '\begin{tabular}{ll} \includegraphics[width=1.9in]{./mugshots/{'.$knots{$localkey}{mugshot}[0].'} & \includegraphics[width=1.9in]{./mugshots/{'.$knots{$localkey}{mugshot}[1].'} \\\\'."\n".'\includegraphics[width=1.9in]{./mugshots/{'.$knots{$localkey}{mugshot}[2].'} & \end{tabular} '."\n\n".'\vspace{-2.8in}'."\n\n";
	
	print BESTIARY '\hspace{0.6in}\includegraphics[width=4.3in,viewport=150 145 430 560]{./strutplots/'.$knots{$localkey}{strutplot}.'} \vfill'."\n\n";

	# Now we print the table text.
	
	print BESTIARY '\begin{center} \begin{tabular}{lllllllll} \toprule'."\n";
	print BESTIARY 'Link & $\PRop$ & $\Rop$ & Filename & Verts & Struts & $\kappa$ range & Kink & Straight \\\\ \midrule'."\n";
	
	print BESTIARY '$'.$texkey.'$ & $'.$knots{$localkey}{ropelength}.'$ & $'.$knots{$localkey}{upperbound}.'$ & \verb!'.$knots{$localkey}{filename}.'! & $'.$knots{$localkey}{edges}.'$ & $'.$knots{$localkey}{struts}.'$ & $['.$knots{$localkey}{klo}.','.$knots{$localkey}{khi}.']$ & $'.$knots{$localkey}{kinks}.'$ & $'.$knots{$localkey}{straight}.'$ \\\\'."\n";
	
	print BESTIARY '\bottomrule'."\n".'\end{tabular} \end{center}'."\n\n";
	
    }
     
     print "\n\n";
     
 }


close BESTIARY;   
    
