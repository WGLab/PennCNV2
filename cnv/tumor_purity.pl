#!/usr/bin/env perl

use warnings;
use strict;
use Carp;
use Getopt::Long;
use Pod::Usage;

our ($verbose, $help, $man);
our ($inputfile, $hmmfile);
our ($output, $snpposfile, $bin, $grid);

GetOptions ('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'output=s'=>\$output, 'snpposfile=s'=>\$snpposfile, 
	'bin=i'=>\$bin, 'grid=i'=>\$grid) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 2 or pod2usage ("Syntax error:");

($inputfile, $hmmfile) = @ARGV;

$bin ||= 100;
$bin > 2 or pod2usage ("Error: the -bin argument must be higher than 2");
$grid ||= 10;
$grid > 2 or pod2usage ("Error in argument: the --grid argument must be higher than 2");

if (defined $output) {
	open (STDOUT, ">$output") or confess "Error: cannot write to output file $output: $!\n";
}

my $siginfo = readSignal ($inputfile);
my $hmm = readHMMFile ($hmmfile);

predictFraction ($siginfo, $hmm);


sub predictFraction {
	my ($siginfo, $hmm) = @_;
	my (@lrr, @baf);
	my (@bin_lrr_maxi);
	for my $key (keys %$siginfo) {
		$key =~ m/^\d+$/ or next;		#igore X, Y, MT, etc.
		for my $i (0 .. @{$siginfo->{$key}}-1) {
			my ($pos, $lrr, $baf) = @{$siginfo->{$key}[$i]};
			$lrr =~ m/^[\d\.\-eE]+$/ or next;
			$baf =~ m/^[\d\.eE]+$/ or next;
			
			
			push @lrr, $lrr;
			push @baf, $baf;
			
			if (@lrr == $bin) {		#reached binning threshold
				push @bin_lrr_maxi, [mean (\@lrr), calculateBAFLikelihood (\@baf, $hmm)];
				@lrr = ();
				@baf = ();
			}
		}
	}
	@bin_lrr_maxi = sort {$a->[0]<=>$b->[0]} @bin_lrr_maxi;
	
	my %histogram;
	for my $i (0 .. @bin_lrr_maxi-1) {		#results do not change much for the all bins, or the first 20% bins with lowerst mean LRR values (deletions)
		$histogram{$bin_lrr_maxi[$i]->[1]}++;
	}
	
	my (@alphacount, $maxalpha, $maxindex);
	for my $i (0 .. $grid) {
		$histogram{$i}||=0;
		print STDERR "$i\t$histogram{$i}\n";
		
		push @alphacount, [$grid, $histogram{$i}];
		if (not defined $maxindex) {
			$maxindex = $i;
			$maxalpha = $histogram{$i};
		}
		if ($histogram{$i} > $maxalpha) {
			$maxindex = $i;
			$maxalpha = $histogram{$i};
		}
	}
	print "$inputfile maxindex = $maxindex, maxcount = $maxalpha tumor_purity(1-alpha) = ", $maxindex/$grid, "\n";
	
	#comment out below (the below calculate the weighted average estimate of the alpha. finding mode is perhaps more reliable than this)
	#my ($sum, $count) = (0, 0);
	#for my $i (0 .. $grid) {
	#	$count += $histogram{$i};
	#	$sum += $histogram{$i}*$i/$grid;
	#}
	#print STDERR "average alpha estimate=", $sum/$count, "\n";
}

sub calculateBAFLikelihood {
	my ($allbaf, $hmm) = @_;
	my @baf = @$allbaf;
	my $h = 0.3;
	my @mu;
	my @sd;
	
	for my $i (0 .. $grid) {
		push @mu, 0.5 * (1-$i/$grid); 		# qw/0.5 0.45 0.4 0.35 0.3 0.25 0.2 0.15 0.1 0.05 0/;
		#push @mu, (1-$i/$grid)/(2-$i/$grid);
	}
	
	my ($sdcn0, $sdcn2) = ($hmm->{'B2_sd'}[0], $hmm->{'B2_sd'}[3]);
	
	for my $i (0 .. $grid) {
		push @sd, $sdcn2 - ($sdcn2-$sdcn0)* $i/$grid;
	}
	
	#print "mu=@mu\n";
	#print "sd=@sd\n";
	my ($lnorm, @ldel);
	for my $i (0 .. @baf-1) {
		$baf[$i] > 0.5 and $baf[$i] = 1-$baf[$i];
		$baf[$i] >=0 and $baf[$i] <=1 or print STDERR "WARNING: skipping invalid BAF values $baf[$i]\n" and next;
		for my $j (0 .. $grid) {
			my $like = (1-$h) * pdf_normal ($baf[$i], 0, $sd[0]) + $h * pdf_normal ($baf[$i], $mu[$j], $sd[$j]);
			if ($like) {
				$ldel[$j] += log ($like);
			} else {
				print STDERR "WARNING: zero likelihood found for baf=$baf[$i] grid=$j mu[j]=$mu[$j] sd[j]=$sd[$j]\n";
				next;
			}
		}
	}
	
	my ($maxi, $maxldel) = (0, $ldel[0]);
	for my $i (1 .. $grid) {
		if ($ldel[$i] >= $maxldel) {
			$maxi = $i;
			$maxldel = $ldel[$i];
		}
	}			
	
	$verbose and print "NOTICE: lnorm=$lnorm bestprediction= $maxi ldel= @ldel\n";
	return ($maxi);
}

#read the HMM file
sub readHMMFile {
	my ($inputfile) = @_;
	my (%hmm, @cell);
	open (HMM, $inputfile) or confess "\nERROR: cannot read from HMM file $hmmfile: $!\n";
	my @line = <HMM>;
	map {s/[\r\n]+$//} @line;
	$line[0] eq 'M=6' or confess "\nERROR: invalid record found in HMM file: <$_> ('M=6' expected)\n";
	$line[1] eq 'N=6' or confess "\nERROR: invalid record found in HMM file: <$_> ('N=6' expected)\n";
	$line[2] eq 'A:' or confess "\nERROR: invalid record found in HMM file: <$_> ('A:' expected)\n";
	$line[9] eq 'B:' or confess "\nERROR: invalid record found in HMM file: <$_> ('B:' expected)\n";
	$line[16] eq 'pi:' or confess "\nERROR: invalid record found in HMM file: <$_> ('pi:' expected)\n";
	$line[18] eq 'B1_mean:' or confess "\nERROR: invalid record found in HMM file: <$_> ('B1_mean:' expected)\n";
	$line[20] eq 'B1_sd:' or confess "\nERROR: invalid record found in HMM file: <$_> ('B1_sd:' expected)\n";
	$line[22] eq 'B1_uf:' or confess "\nERROR: invalid record found in HMM file: <$_> ('B1_uf:' expected)\n";
	$line[24] eq 'B2_mean:' or confess "\nERROR: invalid record found in HMM file: <$_> ('B2_mean:' expected)\n";
	$line[26] eq 'B2_sd:' or confess "\nERROR: invalid record found in HMM file: <$_> ('B2_sd:' expected)\n";
	$line[28] eq 'B2_uf:' or confess "\nERROR: invalid record found in HMM file: <$_> ('B2_uf:' expected)\n";
	if (@line > 30) {
		#print STDERR "NOTICE: HMM model file $inputfile contains parameters for non-polymorphic (NP) probes\n";
		@line == 36 or @line == 38 or confess "\nERROR: invalid number of records found in HMM file: 30 or 36 or 38 lines expected\n";
		$line[30] eq 'B3_mean:' or confess "\nERROR: invalid record found in HMM file: <$_> ('B3_mean:' expected)\n";
		$line[32] eq 'B3_sd:' or confess "\nERROR: invalid record found in HMM file: <$_> ('B3_sd:' expected)\n";
		$line[34] eq 'B3_uf:' or confess "\nERROR: invalid record found in HMM file: <$_> ('B3_uf:' expected)\n";
		if (@line == 38) {
			$line[36] eq 'DIST:' or confess "\nERROR: invalid record found in HMM file: <$_> ('DIST:' expected)\n";
		}
	}

	for my $i (3 .. 8) {
		@cell = split (/\s+/, $line[$i]);
		abs (sum (\@cell) - 1) < 1e-5 or confess "\nERROR: invalid line ${\($i+1)} in HMM file: <$_> (sum of line should be 1)\n";
		push @{$hmm{'A'}}, [@cell];
	}
	
	@cell = split (/\s+/, $line[17]);
	abs (sum (\@cell) - 1) < 1e-5 or confess "\nERROR: invalid line in HMM file: <$line[17]> (sum of line should be 1)\n";
	push @{$hmm{'pi'}}, @cell;
	
	@cell = split (/\s+/, $line[19]);
	@cell == 6 or confess "\nERROR: invalid line in HMM file: <@cell> (6 fields expected)\n";
	push @{$hmm{'B1_mean'}}, @cell;
	
	@cell = split (/\s+/, $line[21]);
	@cell == 6 or confess "\nERROR: invalid line in HMM file: <@cell> (6 fields expected)\n";
	push @{$hmm{'B1_sd'}}, @cell;
	grep {$_>0} @cell == 6 or confess "\nERROR: invalid line in HMM file: <@cell> (all values should be between greater than zero)\n";
	
	@cell = split (/\s+/, $line[23]);
	@cell == 1 or confess "\nERROR: invalid line in HMM file: <@cell> (1 fields expected)\n";
	push @{$hmm{'B1_uf'}}, @cell;
	
	@cell = split (/\s+/, $line[25]);
	@cell == 5 or confess "\nERROR: invalid line in HMM file: <@cell> (5 fields expected)\n";
	push @{$hmm{'B2_mean'}}, @cell;
	
	@cell = split (/\s+/, $line[27]);
	@cell == 5 or confess "\nERROR: invalid line in HMM file: <@cell> (5 fields expected)\n";
	push @{$hmm{'B2_sd'}}, @cell;
	grep {$_>0} @cell == 5 or confess "\nERROR: invalid line in HMM file: <@cell> (all values should be between greater than zero)\n";

	@cell = split (/\s+/, $line[29]);
	@cell == 1 or confess "\nERROR: invalid line in HMM file: <@cell> (1 fields expected)\n";
	push @{$hmm{'B2_uf'}}, @cell;
	
	if (@line > 30) {
		@cell = split (/\s+/, $line[31]);
		@cell == 6 or confess "\nERROR: invalid line in HMM file: <@cell> (6 fields exepcted)\n";
		push @{$hmm{'B3_mean'}}, @cell;

		@cell = split (/\s+/, $line[33]);
		@cell == 6 or confess "\nERROR: invalid line in HMM file: <@cell> (6 fields expected)\n";
		push @{$hmm{'B3_sd'}}, @cell;
		grep {$_>0} @cell == 6 or confess "\nERROR: invalid line in HMM file: <@cell> (all values should be between greater than zero)\n";

		@cell = split (/\s+/, $line[35]);
		@cell == 1 or confess "\nERROR: invalid line in HMM file: <@cell> (1 fields expected)\n";
		push @{$hmm{'B3_uf'}}, @cell;
	}
	
	if (@line == 38) {
		$line[37] =~ m/^\d+$/ or confess "Error: invalid line in HMM file: <$line[37]> (an integer expected for DIST)\n";
		$hmm{dist} = $line[37];
	}
	
	close (HMM);
	return (\%hmm);
}		


sub readSignal {
	my ($signalfile) = @_;
	my (%bin, %siginfo);
	my $countsnp = 0;
	
	open (FH, $signalfile) or die "Error: cannot read from signalfile $signalfile: $!\n";
	print STDERR "NOTICE: Reading signal intensity information from $signalfile ...";
	$_ = <FH>;
	defined ($_) or die "\nERROR: NOTHING is found in signalfile $signalfile. Please check the file before proceeding\n";
	s/[\r\n]+$//;
	my @header = split (/\t/, $_);
	my ($name_index, $chr_index, $pos_index, $lrr_index, $baf_index);
	for my $i (0 .. @header-1) {
		if ($header[$i] eq 'Name' or $header[$i] eq 'SNP ID') {
			$name_index = $i;
		} elsif ($header[$i] eq 'Chr' or $header[$i] eq 'Chromosome') {
			$chr_index = $i;
		} elsif ($header[$i] eq 'Position') {
			$pos_index = $i;
		} elsif ($header[$i] =~ m/Log R Ratio$/ or $header[$i] eq 'LRR') {
			$lrr_index = $i;
		} elsif ($header[$i] =~ m/B Allele Freq$/ or $header[$i] eq 'BAF') {
			$baf_index = $i;
		}
	}

	defined $name_index or die "Error: the signal file $signalfile does not contain the Name column in the header line\n";
	if (not defined $snpposfile) {
		defined $chr_index and defined $pos_index or die "Error: the signal file $signalfile does not contain Chr and Position column, use --snpposfile to supply this information\n";
	}
	defined $name_index and defined $lrr_index and defined $baf_index or die "Error: the signalfile $signalfile does not contain the Name, LRR and BAF columns in the first line\n";
	
	my ($name_chr, $name_pos);
	$snpposfile and ($name_chr, $name_pos) = readSNPPosFile ($snpposfile);
	
	while (<FH>) {
		s/[\r\n]+$//;
		my @field = split (/\t/, $_);
		my ($name, $chr, $pos, $lrr, $baf);
		if (defined $chr_index and defined $pos_index) {
			($name, $chr, $pos, $lrr, $baf) = @field[$name_index, $chr_index, $pos_index, $lrr_index, $baf_index];
		} else {
			($name, $lrr, $baf) = @field[$name_index, $lrr_index, $baf_index];
			($chr, $pos) = ($name_chr->{$name}, $name_pos->{$name});
		}
		defined $name and defined $chr and defined $pos and defined $lrr and defined $baf or next;		#this record does not contain required information
		
		
		
		#$name =~ m/cnv/ and next;			#exclude cnv markers from analysis (for Illumina arrays)
		#$name =~ m/CN/ and next;			#exclude cnv markers from analysis (for Illumina arrays)
		#$name =~ m/^rs/ or $name =~ m/^SNP/ or next;	#allow users to use custom names
		
		
		$countsnp++;
		$lrr =~ m/^[\d\.\-\e]+$/ or next;		#make sure that they are numbers
		$baf =~ m/^[\d\.\e]+$/ or next;
				
		push @{$siginfo{$chr}}, [$pos, $lrr, $baf];
	}
	
	for my $key (keys %siginfo) {		#sort by chromosome position
		@{$siginfo{$key}} = sort {$a->[0] <=> $b->[0]} @{$siginfo{$key}};
	}
	print STDERR "Done with $countsnp SNP markers in ", scalar (keys %siginfo), " chromosomes.\n";
	return (\%siginfo);
}

sub readSNPPosFile {
	my ($snpposfile) = @_;
	my (%name_pos, %name_chr);
	my ($header, $header_seg, $lrr_index, $pos_index, $chr_index, $name_index);
	
	open (SIG, $snpposfile) or die "Error: cannot read from SNPLOC file $snpposfile: $!\n";
	print STDERR "NOTICE: Start reading snpposfile $snpposfile ...";
	$header = <SIG>;
	$header =~ m/(.+)Pos/ or die "error: the header line of $snpposfile does not contain Pos annotation";
	$header_seg = $1;
	$pos_index = ($header_seg =~ tr/\t/\t/+0);
	$header =~ m/(.+)Chr/ or die "error: the header file of $snpposfile does not contain Chr annotation";
	$header_seg = $1;
	$chr_index = ($header_seg =~ tr/\t/\t/+0);
	$header =~ m/(.*)Name/ or die "error: the header file of $snpposfile does not contain Name annotation";
	$header_seg = $1;
	$name_index = ($header_seg =~ tr/\t/\t/+0);


	while (<SIG>) {
		s/[\r\n]+$//;
		my @record = split (/\t/, $_);
		my ($curname, $curchr, $curpos) = @record[$name_index, $chr_index, $pos_index];
		$name_chr{$curname} = $curchr;
		$name_pos{$curname} = $curpos;
	}
	close (SIG);
	print STDERR " Done with location information for ${\(scalar keys %name_chr)} markers\n";
	return (\%name_chr, \%name_pos);
}

#the following subroutine calculates the correlation coefficient
sub cc {
	my ($score1, $score2) = @_;
	my ($ssr12, $ssr11, $ssr22) = (ssr ($score1, $score2), ssr ($score1, $score1), ssr ($score2, $score2));
	$ssr11*$ssr22 or return "NA";
	return $ssr12 / sqrt ($ssr11 * $ssr22);
}

#the following subroutine calculates the ssr score, which is used in cc (correlation coefficient) calculation
sub ssr {
	my @score1 = @{$_[0]};
	my @score2 = @{$_[1]};
	my $mean1 = mean ($_[0]);
	my $mean2 = mean ($_[1]);
	my $product = 0;
	for my $i (0 .. @score1-1) {
		$product += $score1[$i] * $score2[$i];
	}
	return ($product - @score1 * $mean1 * $mean2);
}	

sub sum {
	my ($score) = @_;
	@$score or confess "\nERROR: NO VALUES for calculating sum\n";
	my $sum;
	for (@$score) {
		$sum += $_;
	}
	return $sum;
}

sub mean {
	my ($score) = @_;
	@$score or confess "\nERROR: NO VALUES for calculating mean\n";
	my $sum;
	for (@$score) {
		$sum += $_;
	}
	return $sum/@$score;
}

sub sd {
	my ($score) = @_;
	@$score > 1 or confess "\nERROR: NO sufficient VALUES for calculating SD\n";
	my $mean = mean ($score);
	my $sum;
	for my $i (0 .. @$score-1) {
		$sum += ($score->[$i]-$mean)*($score->[$i]-$mean);
	}
	$sum /= (@$score-1);
	return sqrt ($sum);
}

sub median {
	my ($score) = @_;
	@$score or confess "\nERROR: NO VALUES for calculating median\n";
	my @newscore = sort {$a<=>$b} @$score;
	if (@newscore % 2 == 0) {
		return ($newscore[@newscore/2-1]+$newscore[@newscore/2])/2;
	} else {
		return $newscore[@newscore/2];
	}
}

sub pdf_normal {
	my ($x, $mu, $sd) = @_;
	return exp(-($x-$mu)*($x-$mu)/2/$sd/$sd)/$sd/sqrt(2*3.1415927);
}
	

=head1 SYNOPSIS

 tumor_cnv.pl [arguments] <input-signal-file> <HMM-file>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	    --snpposfile <file>		a file with chr/position information for markers
 	    --bin <int>			the BIN for grouping SNPs together (default: 100)
 	    --grid <int>		the GRID for precision of alpha estimate (default: 10)
 	

 Function: identify copy number for a mixed population of cells with two components

 Example: tumor_cnv.pl signal.txt hhall.hmm

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=item B<--output>

specify the output file name (default is to print to STDOUT)

=item B<--snpposfile>

a file with Chr and Position information for each marker (if the signalfile does
not contain this information)

=item B<--bin>

the BIN for grouping SNPs together, with a default value of 100 which works well
for typical high-density SNP arrays

=item B<--grid>

the GRID for precision of alpha estimate, with a default value of 10 but a higher
value such as 50 can be specified for high-density SNP arrays

=back

=head1 DESCRIPTION

This program is designed to identify copy number for a mixed population of cells 
with two components. Unlike typical CNV calling programs that aim to give a copy 
number estimates as integers, this program is suitable for samples that contain 
two sub-clones, that is, two populations of cells are mixed together. The copy 
number will not be an integer, but a fractional measure determined by copy 
number of each subclone as well as the relative fraction of each subclone in the 
sample.

There are two main usages of the program:

1. Given the known lack of copy number changes, it can be used to identify 
"fractional LOH" patterns of signal intensity. That is, one subclone is normal, 
yet the subclone has a LOH, resulting in highly characteristic BAF patterns 
across a specific genomic region. This is useful, for example, in calculating 
fractional LOH in 9p in patients with polycythemia vera, a disease known to be 
related to 9p LOH. When using this functionality, it is necessary to always 
specify a specific genomic region.

2. Assuming that a sample is composed of normal subclone and another subclone 
with copy number changes in many genomic regions (such as a tumor tissue with 
stromal contamination), this program can scan whole genome or a specific region 
and try to guess the extent of stromal contamination.

For questions, comments or bug reports, please contact me at 
kai@openbioinformatics.org.

=cut
