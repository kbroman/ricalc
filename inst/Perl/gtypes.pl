#!/usr/local/bin/perl
######################################################################
# gtypes.pl
#
# the goal here is to determine all possible parental types
# for the process of producing RILs by sib mating, and to collapse 
# them into equivalence classes by various symmetries
#
######################################################################

(@ARGV > 2) or die("Give no. strains, chr type, and no. loci\n");

$n_strains = $ARGV[0];
$chrtype = $ARGV[1];
$n_loci = $ARGV[2];

if($n_strains!=2 and $n_strains!=4) {
    die("no. strains should be 2 or 4\n");
}
if($chrtype ne "A" and $chrtype ne "X") {
    die("chr type should be A/X\n");
}
if($n_loci!=2 and $n_loci!=3) {
    die("No. loci should be 2 or 3\n");
}


@alleles = ("A","B","C","D");
if($chrtype eq "X" and $n_strains==4) {
    @alleles = @alleles[0..2];
}
else {
    @alleles = @alleles[0..($n_strains-1)];
}
$na = @alleles;
$nam1 = $na - 1;

@haps = all_hap();
$n_haps = @haps;
print(" No. haplotypes  = $n_haps\n");

@inds = all_ind();
$n_inds = @inds;
print(" No. individuals = $n_inds\n");

print(" -Creating all pairs\n");
if($chrtype eq "X") {
    $s = 0;
    foreach $i (0..($n_inds-1)) {
	foreach $j (0..($n_haps-1)) {
	    $pair[$s] = $inds[$i] . " x " . $haps[$j];
	    $s++;
	}
    }
}
else {
    $s = 0;
    foreach $i (0..($n_inds-1)) {
	foreach $j ($i..($n_inds-1)) {
	    $pair[$s] = $inds[$i] . " x " . $inds[$j];
	    $s++;
	}
    }
}

$n_pair = @pair;
print("No. pairs = $n_pair\n");

print(" -Reducing pairs by symmetry\n");
$n_lookup=0;
foreach $i (0..($n_pair-1)) {
    if($i == sprintf("%d", $i/1000)*1000) {
	printf("%7d of %7d / %d\n", $i, $n_pair, $n_lookup);
    }
    if($lookup{$pair[$i]} eq "") {
	$n_lookup++;
	@alt = alternates($pair[$i]);
	foreach $alt (@alt) {
	    $lookup{$alt} = $pair[$i];
	}
    }
}

print("Number by symmetry: $n_lookup\n");

print(" -Writing output file.\n");
$ofile = "out" . $n_strains . $chrtype . $n_loci . ".csv";
open(OUT, ">$ofile") or die("Cannot write to $ofile");
print OUT ("all,proto\n");
foreach $pair (@pair) {
    print OUT ("$pair,$lookup{$pair}\n");
}
close(OUT);

######################################################################
######################################################################


######################################################################
# subroutine for exchanging alleles
# 
# first argument in input should be a pattern like "AB"
######################################################################
sub exchangeg 
{
    @exch_out = @_;
    $exch_pat = @exch_out[0];
    @exch_pat = split(//, $exch_pat);
    foreach $ei (1..(@exch_out-1)) {
	if($exch_pat eq "AB") {
	    $exch_out[$ei] =~ s/A/Z/g;
	    $exch_out[$ei] =~ s/B/A/g;
	    $exch_out[$ei] =~ s/Z/B/g;
	}
	elsif($exch_pat eq "CD") {
	    $exch_out[$ei] =~ s/C/Z/g;
	    $exch_out[$ei] =~ s/D/C/g;
	    $exch_out[$ei] =~ s/Z/D/g;
	}
	elsif($exch_pat eq "AC") {
	    $exch_out[$ei] =~ s/C/Z/g;
	    $exch_out[$ei] =~ s/A/C/g;
	    $exch_out[$ei] =~ s/Z/A/g;
	}
	elsif($exch_pat eq "BD") {
	    $exch_out[$ei] =~ s/B/Z/g;
	    $exch_out[$ei] =~ s/D/B/g;
	    $exch_out[$ei] =~ s/Z/D/g;
	}
	else {
	    die("Need to program $exch_pat\n");
	}
	    
    }

    $n = @exch_out;
    @exch_out[(1..(@exch_out-1))];
}

	
######################################################################
# subroutine for reversing the locus order of a bunch of patterns
######################################################################
sub reverseg 
{
    @rev_out = @_;
    foreach $ri (0..(@_-1)) {
	@rev_v = split("[\\|x]", $rev_out[$ri]);
	foreach $rj (0..(@rev_v-1)) {
	    $rev_v[$rj] =~ s/\s+//g;
	    $rev_v[$rj] = join("",reverse(split(//, $rev_v[$rj])));
	}
	$nv = @rev_v;
	if($nv == 1) { 
	    $rev_out[$ri] = $rev_v[0]; 
	}
	elsif($nv == 2) { 
	    $rev_out[$ri] = $rev_v[0] . "|" . $rev_v[1]; 
	}
	elsif($nv == 3) { 
	    $rev_out[$ri] = $rev_v[0] . "|" . $rev_v[1] . " x " . $rev_v[2]; 
	}
	else { 
	    $rev_out[$ri] = $rev_v[0] . "|" . $rev_v[1] . " x " . 
		$rev_v[2] . "|" . $rev_v[3]; 
	}
    }
    @rev_out;
}


######################################################################
# adjust_order
#
# Adjust the order of haplotypes / individuals in a pattern
######################################################################
sub adjust_order
{
    @ao_out = @_;

    foreach $aoi (0..(@ao_out-1)) {
	@ao_v = split("[\\|x]", $ao_out[$aoi]);
	foreach $aoj (0..(@ao_v-1)) {
	    $ao_v[$aoj] =~ s/\s+//g;
	}
	$nv = @ao_v;
	if($nv > 1) {
	    @ao_v[0..1] = sort(@ao_v[0..1]);
	}
	if($nv > 3) {
	    @ao_v[2..3] = sort(@ao_v[2..3]);
	}

	if($nv == 1) { 
	    $ao_out[$aoi] = $ao_v[0]; 
	}
	elsif($nv == 2) { 
	    $ao_out[$aoi] = $ao_v[0] . "|" . $ao_v[1]; 
	}
	elsif($nv == 3) { 
	    $ao_out[$aoi] = $ao_v[0] . "|" . $ao_v[1] . " x " . $ao_v[2]; 
	}
	else { 
	    ($a,$b) = sort($ao_v[0] . "|" . $ao_v[1],
			   $ao_v[2] . "|" . $ao_v[3]);
	    $ao_out[$aoi] = $a . " x " . $b;
	}
    }
    @ao_out;
}

######################################################################
# alternates
#
# Uses symmetry to find all patterns equivalent to the input one.
######################################################################
sub alternates 
{
    $ap_pat = $_[0];

    if($n_strains==2 and $chrtype eq "X") {
	@ap_out = adjust_order($ap_pat,reverseg($ap_pat));
    }
    elsif($n_strains==2 or ($n_strains==4 and $chrtype eq "X")) {
	$ap_pat2 = exchangeg("AB",$ap_pat);
	@ap_out = adjust_order($ap_pat,$ap_pat2,
			       reverseg($ap_pat),reverseg($ap_pat2));
    }
    else {
	$ap_pat2 = exchangeg("AB", $ap_pat);
	$ap_pat3 = exchangeg("CD", $ap_pat);
	$ap_pat4 = exchangeg("CD", $ap_pat2);
	$ap_pat5 = exchangeg("BD", exchangeg("AC", $ap_pat));
	$ap_pat6 = exchangeg("AB", $ap_pat5);
	$ap_pat7 = exchangeg("CD", $ap_pat6);
	$ap_pat8 = exchangeg("CD", $ap_pat5);
	@ap_out = adjust_order($ap_pat,$ap_pat2,$ap_pat3,$ap_pat4,
			       $ap_pat5,$ap_pat6,$ap_pat7,$ap_pat8,
			       reverseg($ap_pat),reverseg($ap_pat2),
			       reverseg($ap_pat3), reverseg($ap_pat4),
			       reverseg($ap_pat5),reverseg($ap_pat6),
			       reverseg($ap_pat7),reverseg($ap_pat8));
    }
    sort_unique(@ap_out);
}

######################################################################
# sort_unique
#
# find unique values and sort them
######################################################################

sub sort_unique 
{
    %su_hash = ();

    foreach $sui (@_) {
	$su_hash{$sui} = 1;
    }

    sort keys %su_hash;
}

######################################################################
# all_hap
#
# all haplotypes
######################################################################

sub all_hap 
{
    if($n_loci==3) {
	$ahs=0;
	foreach $ahi (0..$nam1) {
	    foreach $ahj (0..$nam1) {
		foreach $ahk (0..$nam1) {
		    $haps[$ahs] = $alleles[$ahi] . $alleles[$ahj] . $alleles[$ahk];
		    $ahs++;
		}
	    }
	}
    }
    else {
	$ahs=0;
	foreach $ahi (0..$nam1) {
	    foreach $ahj (0..$nam1) {
		$haps[$ahs] = $alleles[$ahi] . $alleles[$ahj];
		$ahs++;
	    }
	}
    }
    @haps;
}

######################################################################
# all_ind
# all individuals
######################################################################
sub all_ind 
{
    $nh = @haps;
    $ais = 0;
    foreach $aii (0..($nh-1)) {
	foreach $aij ($aii..($nh-1)) {
	    $inds[$ais] = $haps[$aii] . "|" . $haps[$aij];
	    $ais++;
	}
    }
    @inds;
}

