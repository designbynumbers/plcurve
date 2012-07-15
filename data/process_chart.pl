open(CHART,"knotchart-prime.txt") or die("Couldn't open knotchart\n");
@lines = <CHART>;
close(CHART);

foreach $line (@lines) {

    $line =~ /(.+)\s(.+)/;
    
    $knot = $1;
    $sym = $2;

    $sym =~ s/reversible/Reversible/;
    $sym =~ s/none/NoSymmetry/;
    $sym =~ s/full/FullSymmetry/;
    $sym =~ s/amphichiral/PlusAmphichiral/;
    $sym =~ s/\s//g;

    $knot =~ s/m//g;
    $knot =~ s/r//g;
    $knot =~ s/\s//g;

    $khash{$knot} = $sym;

}

foreach $key (sort {

    knot($a) <=> knot($b) or

	ind($a) <=> ind($b)

	      } keys %khash) {

    print "\{\"$key\",$khash{$key},".knot($key).",".ind($key)."\},\n";

}


sub knot {

    my $string = shift;
    @parts = split("_",$string);
    return $parts[0];

}

sub ind {

    my $string = shift;
    @parts = split("_",$string);
    return $parts[1];

}


