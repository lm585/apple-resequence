#!/usr/bin/perl -w
# Used to check duplicated reads in sweetpotato reads.
# Edit for WM_D1_test 
# Edit 2012-10-09
use strict;
use Getopt::Long; 


my $ptag = 'out'; 
my $help = ''; 
my $is_ogz = ''; 
my $no_out = ''; 
GetOptions(
'opre=s' => \$ptag, 
'help' => \$help, 
'ogz' => \$is_ogz,
'noout' => \$no_out, 
); 

( !@ARGV or $help ) and die "perl $0 -ogz -opre out_prefix input1.fastq input2.fastq\n"; 
(-f $ARGV[0] and -f $ARGV[1]) or die "perl $0 -ogz -opre out_prefix input1.fastq input2.fastq\n"; 

my ($rf1, $rf2) = @ARGV; 

if ($rf1 =~ m/\.gz$/) {
	open (R1,'-|', "gzip -cd $rf1") or die; 
}else{
	open (R1, '<', "$rf1") or die; 
}
if ($rf2 =~ m/\.gz$/) {
	open (R2,'-|', "gzip -cd $rf2") or die; 
}else{
	open (R2, '<', "$rf2") or die; 
}

unless ($no_out) {
	if ($is_ogz) {
		open (O1,'|-', "gzip -c >${ptag}_R1.ndupB.gz") or die;
		open (O2,'|-', "gzip -c >${ptag}_R2.ndupB.gz") or die; 
		open (D1,'|-', "gzip -c >${ptag}_R1.dropB.gz") or die; 
		open (D2,'|-', "gzip -c >${ptag}_R2.dropB.gz") or die;
		warn "[Stat]Reading $ptag files: $rf1 $rf2\n"; 
		warn "[Stat]Create output files: ${ptag}_R1.ndupB.gz ${ptag}_R2.ndupB.gz ${ptag}_R1.dropB.gz ${ptag}_R2.dropB.gz\n"; 
	}else{
		open (O1,'>', "${ptag}_R1.ndupB") or die;
		open (O2,'>', "${ptag}_R2.ndupB") or die;
		open (D1,'>', "${ptag}_R1.dropB") or die; 
		open (D2,'>', "${ptag}_R2.dropB") or die; 
		warn "[Stat]Reading $ptag files: $rf1 $rf2\n"; 
		warn "[Stat]Create output files: ${ptag}_R1.ndupB ${ptag}_R2.ndupB ${ptag}_R1.dropB ${ptag}_R2.dropB\n"; 
	}
}

my (%have1, %have2); 
my %have_both; 
my $total = 0; 
my $uniq = 0; 
my $mult = 0; 
my $uniq1 = 0; my $uniq2 = 0; 
my $mult1 = 0; my $mult2 = 0; 
while (!eof(R1)) {
	$. % 10000000 == 0 and warn "[Stat]$. line in $ptag files. $mult duplicated read pairs now.\n"; 
	my ($k1, $s1, $q1) = &readArecord(\*R1);
	my ($k2, $s2, $q2) = &readArecord(\*R2);
	$total ++; 
	my $s_both = "$s1\t$s2"; 
	if (defined $have_both{$s_both}) {
		if ($have_both{$s_both} == 1) {
			$mult += 2; 
		}else{
			$mult ++; 
		}
		$have_both{$s_both} ++; 
		if ($have1{$s1} == 1) {
			$mult1 += 2; 
		}else{
			$mult1 ++; 
		}
		$have1{$s1} ++; 
		if ($have2{$s2} == 1) {
			$mult2 += 2; 
		}else{
			$mult2 ++; 
		}
		$have2{$s2} ++; 
	unless ($no_out) {
			print D1 "$k1 [BDup]\n$s1\n+\n$q1\n"; 
			print D2 "$k2 [BDup]\n$s2\n+\n$q2\n"; 
	}
	}else{
		$uniq ++; 
		$have_both{$s_both} = 1; 
	unless ($no_out) {
			print O1 "$k1\n$s1\n+\n$q1\n"; 
			print O2 "$k2\n$s2\n+\n$q2\n"; 
	}
		if (defined $have1{$s1}) {
			if ($have1{$s1} == 1) {
				$mult1 += 2; 
			}else{
				$mult1++; 
			}
			$have1{$s1}++; 
		}else{
			$uniq1++; 
			$have1{$s1} = 1; 
		}
		if (defined $have2{$s2}) {
			if ($have2{$s2} == 1) {
				$mult2 += 2; 
			}else{
				$mult2 ++; 
			}
			$have2{$s2}++; 
		}else{
			$uniq2++; 
			$have2{$s2} = 1; 
		}
	}
}

unless ($no_out) {
	close D1; 
	close D2; 
	close O1;
	close O2;
}
close R1;
close R2;
my $mult_pat_b = $uniq+$mult-$total; 
my $uniq_pat_b = $uniq-$mult_pat_b; 
my $drop_rd_b = $mult-$mult_pat_b; 
my $mult_pat_bR = &rate($mult_pat_b, $total); 
my $uniq_pat_bR = &rate($uniq_pat_b, $total); 
my $drop_rd_bR = &rate($drop_rd_b, $total); 

my $mult_pat_1 = $uniq1+$mult1-$total; 
my $uniq_pat_1 = $uniq1-$mult_pat_1; 
my $drop_rd_1 = $mult1-$mult_pat_1; 
my $mult_pat_1R = &rate($mult_pat_1, $total); 
my $uniq_pat_1R = &rate($uniq_pat_1, $total); 
my $drop_rd_1R = &rate($drop_rd_1, $total); 

my $mult_pat_2 = $uniq2+$mult2-$total; 
my $uniq_pat_2 = $uniq2-$mult_pat_2; 
my $drop_rd_2 = $mult2-$mult_pat_2; 
my $mult_pat_2R = &rate($mult_pat_2, $total); 
my $uniq_pat_2R = &rate($uniq_pat_2, $total); 
my $drop_rd_2R = &rate($drop_rd_2, $total); 

warn "[Stat]Finish $ptag files.\n"; 
warn "[Record] There are $total read pairs in total [$ptag].\n"; 
warn "[Record] There are $uniq_pat_b ($uniq_pat_bR\%) unique patterns, $mult_pat_b ($mult_pat_bR\%) multiple patterns, and $drop_rd_b ($drop_rd_bR\%) reads dropped in both.\n"; 
warn "[Record] There are $uniq_pat_1 ($uniq_pat_1R\%) unique patterns, $mult_pat_1 ($mult_pat_1R\%) multiple patterns, and $drop_rd_1 ($drop_rd_1R\%) reads dropped in 1st file.\n";
warn "[Record] There are $uniq_pat_2 ($uniq_pat_2R\%) unique patterns, $mult_pat_2 ($mult_pat_2R\%) multiple patterns, and $drop_rd_2 ($drop_rd_2R\%) reads dropped in 2nd file.\n";

sub readArecord {
		my ($fh) = @_;
		my $tk = <$fh>;
		my $ts = <$fh>;
		<$fh>;
		my $tq = <$fh>;
		$tk =~ s/\s+$//; 
		$ts =~ s/\s+$//;
		$tq =~ s/\s+$//;
		return ($tk, $ts, $tq);
}

sub rate {
	my ($n1, $n2) = @_; 
	if ($n2 == 0) {
		warn "Denominator is zero.\n"; 
		return 'NA'; 
	}else{
		my $r = int($n1/$n2*10000+0.5)/100; 
		return $r; 
	}
}
