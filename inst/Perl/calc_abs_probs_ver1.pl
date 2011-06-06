#!/usr/local/bin/perl
######################################################################
# calc_abs_probs.pl
#
# the goal here is to calculate the absorbtion probabilities
# for 4-way RILs by sibling mating, for the autosome or X chromosome,
# and for a particular recombination fraction and coincidence.
#
######################################################################

(@ARGV > 1) or die("Give chr type and file with  rec frac, and coincidences\n");

$chrtype = $ARGV[0];
$rfile = $ARGV[1];
$n_loci = 3;
$n_strains = 4;
$tol = 1e-14;
$maxit = 1000;

$| = 1;

print(" -Reading recombination fractions\n");
open(IN, $rfile) or die("Cannot read from $rfile");
$line = <IN>;
while($line = <IN>) {
    chomp($line);
    ($r, $coi) = split(/,/, $line);
    push(@allr, $r);
    push(@allcoi, $coi);
}
close(IN);

if($chrtype ne "A" and $chrtype ne "X") {
    die("chr type should be A or X\n");
}

if($chrtype eq "X") {
    @alleles = ("A","B","C");
}
else {
    @alleles = ("A","B","C","D");
}

$na = @alleles;
$nam1 = $na - 1;

print(" -Reading lookup table\n");
$file = "lookup4" . $chrtype . "3.csv";
unless(-e $file) {
    $rezip = 1;
    system("gunzip $file");
}
open(IN, $file) or die("Cannot read from $file");
$line = <IN>;
while($line = <IN>) {
    chomp($line);
    ($all,$proto) = split(/,/, $line);
    $lookup{$all} = $proto;
    $states{$proto} = 1;
}
@ustates = sort keys %states;

close(IN);
if($rezip) {
    system("gzip $file");
}

print(" -Finding starting state and absorbing states\n");
if($chrtype eq "X") {
    $start = "AAA|BBB x CCC";
    @absorb = ("AAA", "AAB", "ABA","AAC","ACC","ACA","CAC","ABC","ACB","CCC");
    foreach $absorb (@absorb) {
	$fabsorb{$absorb} = $absorb . "|" . $absorb . " x " . $absorb;
    }
}
else {
    $start = "AAA|BBB x CCC|DDD";
    @absorb = ("AAA", "AAB", "ABA","AAC","ACA","ABC","ACB");
    foreach $absorb (@absorb) {
	$fabsorb{$absorb} = $absorb . "|" . $absorb . " x " . 
	    $absorb . "|" . $absorb;
    }
}
	
print(" -Determining haplotypes and individuals\n");
all_hap();
all_ind();

@temp = split(/\./, $rfile);
$ofile = $temp[0] . "_out." . $temp[1];
open(OUT, ">$ofile") or die("Cannot write to $ofile");

# disable buffering
$old_fh = select(OUT);
$| = 1;
select($old_fh);

print OUT ("rf,coi");
foreach $state (@absorb) {
    print OUT (",$state");
}
print OUT ("\n");

$nr = @allr;
foreach $k (0..(@allr-1)) {
    $rf = $allr[$k];
    $coi = $allcoi[$k];

    printf("rf = %.8f     coi = %.8f     (%d of %d)\n", $rf, $coi, $k+1, $nr);
    print(" -Calculating meiosis matrix\n");
    get_mtm();

    print(" -Calculating absorption probabilities\n");
    %curp = ();
    $curp{$start} = 1;
    foreach $iter (1..$maxit) {
	%nextp = ();


	@curstate = ();
	foreach $state (keys %curp) {
	    if($curp{$state} > $tol/100) { 
		push(@curstate, $state);
	    }
	}
	$n = @curstate;

	printf(" ---iteration %4d  %7d states\n", $iter, $n);

	foreach $state (@curstate) {
	    ($mom, $dad) = split(/ x /, $state);
	    
	    @kids = @prob = ();
	    if($chrtype eq "X") {
		foreach $i (keys %{$mtm{$mom}}) {
		    foreach $j (keys %{$mtm{$mom}}) {
			push(@kids, $i . "|" . $dad .  " x " . $j);
			push(@prob, $mtm{$mom}{$i} * $mtm{$mom}{$j});
		    }
		}
	    }
	    else {
		foreach $i1 (keys %{$mtm{$mom}}) {
		    foreach $i2 (keys %{$mtm{$mom}}) {
			foreach $i3 (keys %{$mtm{$dad}}) {
			    foreach $i4 (keys %{$mtm{$dad}}) {
				push(@kids, $i1 . "|" . $i3 .  " x " . $i2 . "|" . $i4);
				push(@prob, $mtm{$mom}{$i1} * $mtm{$mom}{$i2} *
				     $mtm{$dad}{$i3} * $mtm{$dad}{$i4});
			    }
			}
		    }
		}
	    }
	
	 
	    @kids = adjust_order(@kids);
	    foreach $i (0..(@kids-1)) {
		$kids[$i] = $lookup{$kids[$i]};
	    }

	    foreach $i (0..(@kids-1)) {
		$nextp{$kids[$i]} += ($curp{$state} * $prob[$i]);
	    }

	} # end loop over states
	
	$flag = 0;
	foreach $state (keys %nextp) {
	    if(abs($nextp{$state}-$curp{$state}) > $tol) {
		$flag = 1;
		last;
	    }
	}
	if(!$flag) { last; }
	%curp = %nextp;

    } # end loop over iterations

    print OUT ("$rf,$coi");
    foreach $state (@absorb) {
	printf OUT (",%.15f",$nextp{$fabsorb{$state}});
    }
    print OUT ("\n");

} # end loop over recombination fractions

close(OUT);
    


######################################################################
# get_mtm
#
# calculate meiosis matrix
######################################################################
sub get_mtm
{
    @p = ((1 - 2*$rf + $coi*$rf*$rf)/2, $rf*(1-$coi*$rf)/2, $coi*$rf*$rf/2);
    @p = @p[(0,0,1,1,1,1,2,2)];

    %mtm = ();

    foreach $ind (@inds) {
	@v = split(/\|/, $ind);
	@g1 = split(//, $v[0]);
	@g2 = split(//, $v[1]);
	
	@egg = ($v[0], $v[1], 
		join("", (@g1[(0,1)],$g2[2])),
		join("", (@g2[(0,1)],$g1[2])),
		join("", ($g1[0],@g2[(1,2)])),
		join("", ($g2[0],@g1[(1,2)])),
		join("", ($g1[0],$g2[1],$g1[2])),
		join("", ($g2[0],$g1[1],$g2[2])));
	
	foreach $i (0..(@egg-1)) {
	    $mtm{$ind}{$egg[$i]} += $p[$i];
	}
    }
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

