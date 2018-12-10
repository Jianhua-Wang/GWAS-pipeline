#!/usr/bin/perl
#
# [ vcf2impute_gen.pl ]
#
# This script takes a VCF file and converts it into a .gen file in IMPUTE format.
#
# *** Copyright Bryan Howie, 2013 ***
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

use strict;
use warnings;

sub InitializeSNPAllelesTable();
sub InitializeGtypeProbsLookupTable();
sub PrintGenFile($$);
sub ReadSampleFile($\%);
sub ParseCommandLine();

ParseCommandLine();

# command line variables
my %arg;
my $vcf_file = $arg{-vcf};                    # VCF file to process; can be phased or unphased
my $gen_file = $arg{-gen};                    # name of new file to print; suffix and gzipping are auto
my $chrom = $arg{-chr};                       # chromosome to include in output files
my $start_pos = $arg{-start};                 # first position to include in output files
my $end_pos = $arg{-end};                     # last position to include in output files
my $samp_include_file = $arg{-samp_include};  # file containing a list of samples to include in gen file
my $samp_exclude_file = $arg{-samp_exclude};  # file containing a list of samples to exclude from gen file
my $is_samp_snptest = $arg{-samp_snptest};    # boolean: should we print the sample file in SNPTEST format?
my $is_special_gt = $arg{-special_gt};        # boolean: allow 'special' genotype coding in VCF?
my $is_gt_like = $arg{-gt_like};              # boolean: print genotype likelihoods rather than hard calls?
my $no_mono = $arg{-no_mono};                 # boolean: omit monomorphic sites from processed files?
my $snps_only = $arg{-snps_only};             # boolean: omit any sites that are not biallelic SNPs?
my $indels_only = $arg{-indels_only};         # boolean: omit any sites that are not biallelic INDELs?
my $svs_only = $arg{-svs_only};               # boolean: omit any sites that are not biallelic SVs?

# global variables
my $line;
my (%snp_alleles, %unphased_gt_to_probs, %phased_gt_to_probs, %special_gt_to_probs, %included_samples, %excluded_samples);

# file paths and temporary file names


MAIN: {

  #
  InitializeSNPAllelesTable();

  #
  InitializeGtypeProbsLookupTable();

  #
  if ($samp_include_file ne "") {
    ReadSampleFile($samp_include_file, %included_samples);
  }

  #
  if ($samp_exclude_file ne "") {
    ReadSampleFile($samp_exclude_file, %excluded_samples);
  }

  #
  PrintGenFile($vcf_file, $gen_file);

}



###################################################################
#                            FUNCTIONS                            #
###################################################################

#
sub InitializeSNPAllelesTable() {
  $snp_alleles{'A'} = 1;
  $snp_alleles{'a'} = 1;
  $snp_alleles{'C'} = 1;
  $snp_alleles{'c'} = 1;
  $snp_alleles{'G'} = 1;
  $snp_alleles{'g'} = 1;
  $snp_alleles{'T'} = 1;
  $snp_alleles{'t'} = 1;
}


#
sub InitializeGtypeProbsLookupTable() {
  $unphased_gt_to_probs{'.'} = '0.33 0.33 0.33'; # not sure about this convention
  $unphased_gt_to_probs{'./.'} = '0.33 0.33 0.33';
  $unphased_gt_to_probs{'0/.'} = '0.33 0.33 0.33'; # this can happen in Complete Genomics data
  $unphased_gt_to_probs{'1/.'} = '0.33 0.33 0.33'; # this can happen in Complete Genomics data
  $unphased_gt_to_probs{'0/2'} = '0.33 0.33 0.33'; # ignore since .gen format supports only biallelic variants
  $unphased_gt_to_probs{'1/2'} = '0.33 0.33 0.33'; # ignore since .gen format supports only biallelic variants
  $unphased_gt_to_probs{'2/2'} = '0.33 0.33 0.33'; # ignore since .gen format supports only biallelic variants
  $unphased_gt_to_probs{'0/3'} = '0.33 0.33 0.33'; # ignore since .gen format supports only biallelic variants
  $unphased_gt_to_probs{'1/3'} = '0.33 0.33 0.33'; # ignore since .gen format supports only biallelic variants
  $unphased_gt_to_probs{'2/3'} = '0.33 0.33 0.33'; # ignore since .gen format supports only biallelic variants
  $unphased_gt_to_probs{'3/3'} = '0.33 0.33 0.33'; # ignore since .gen format supports only biallelic variants
  $unphased_gt_to_probs{'0/0'} = '1 0 0';
  $unphased_gt_to_probs{'0/1'} = '0 1 0';
  $unphased_gt_to_probs{'1/0'} = '0 1 0';
  $unphased_gt_to_probs{'1/1'} = '0 0 1';

  $phased_gt_to_probs{'.|.'} = '0.33 0.33 0.33';
  $phased_gt_to_probs{'0|0'} = '1 0 0';
  $phased_gt_to_probs{'0|1'} = '0 1 0';
  $phased_gt_to_probs{'1|0'} = '0 1 0';
  $phased_gt_to_probs{'1|1'} = '0 0 1';

  $special_gt_to_probs{'.'} = '1 0 0';
  $special_gt_to_probs{'./.'} = '1 0 0';
}


#
sub ReadSampleFile($\%) {
  my ($samp_file, $samp_tbl) = @_;
  my @line_info;

  open(SAMP, $samp_file) || die $!;

  while (defined ($line = <SAMP>)) {
    @line_info = split /\s+/, $line;
    $$samp_tbl{$line_info[0]} = 1;
  }

  close(SAMP);
}


#
sub PrintGenFile($$) {
  my ($infile, $outfile) = @_;

  # if the supplied output file name (.gen format) already includes a '.gz'
  # suffix, remove it here
  my $is_out_gz = ($outfile =~ m/.gz$/);
  if ($outfile =~ m/.gz$/) {
    $outfile = substr($outfile, 0, length($outfile)-3);
  }

  my $samp_file = $outfile.".samples";
  $outfile .= ".gz";
  my @line_info;

  # open input VCF file, allowing for possible gzipping
  my $is_in_gz = ($infile =~ m/.gz$/);
  open(VCF, ($is_in_gz ? "gzip -dc $infile |" : $infile)) || die "Can't open $infile: $!\n";

  # create a handle for writing a new gzipped .gen file
  open(GEN, " | gzip -c > $outfile") || die "Can't open $outfile: $!\n";

  # find last header line and use it to print a sample list
  my $pos_tag_idx = -1;
  my $id_tag_idx = -1;
  my $ref_tag_idx = -1;
  my $alt_tag_idx = -1;
  my $filter_tag_idx = -1;
  my $format_tag_idx = -1;
  my @sample_keep_inds = ();

  while (defined ($line = <VCF>)) {
    chomp $line;
    @line_info = split /\t/, $line;
    if (scalar(@line_info) > 0 && $line_info[0] eq "#CHROM") {
      for my $i (1..$#line_info) {
	if ($line_info[$i] eq "POS") {
	  $pos_tag_idx = $i;
	  next;
	}

	if ($line_info[$i] eq "ID") {
	  $id_tag_idx = $i;
	  next;
	}

	if ($line_info[$i] eq "REF") {
	  $ref_tag_idx = $i;
	  next;
	}

	if ($line_info[$i] eq "ALT") {
	  $alt_tag_idx = $i;
	  next;
	}

	if ($line_info[$i] eq "FILTER") {
	  $filter_tag_idx = $i;
	  next;
	}

	if ($line_info[$i] eq "FORMAT") {
	  $format_tag_idx = $i;
	  last;
	}
      }

      # no POS tag found; something must be wrong
      if ($pos_tag_idx == -1) {
	die "\nERROR: No POS tag found.\n\n";
      }

      # no ID tag found; something must be wrong
      if ($id_tag_idx == -1) {
	die "\nERROR: No ID tag found.\n\n";
      }

      # no REF tag found; something must be wrong
      if ($ref_tag_idx == -1) {
	die "\nERROR: No REF tag found.\n\n";
      }

      # no ALT tag found; something must be wrong
      if ($alt_tag_idx == -1) {
	die "\nERROR: No ALT tag found.\n\n";
      }

      # no FILTER tag found; something must be wrong
      if ($filter_tag_idx == -1) {
	die "\nERROR: No FILTER tag found.\n\n";
      }

      # no FORMAT tag found; something must be wrong
      if ($format_tag_idx == -1) {
	die "\nERROR: No FORMAT tag found.\n\n";
      }

      # if we reach this point, we must have successfully found the FORMAT tag;
      # all subsequent entries should be sample IDs
      open(SAMP, "> $samp_file") || die $!;
      if ($is_samp_snptest) { # print SNPTEST header information if requested
	print SAMP "ID_1 ID_2 missing\n";
	print SAMP "0 0 0\n"
      }
      for my $j (($format_tag_idx+1)..$#line_info) {
	if ($samp_include_file ne "" && !exists($included_samples{$line_info[$j]})) {
	  push @sample_keep_inds, 0;
	  next;
	}
	if ($samp_exclude_file ne "" &&  exists($excluded_samples{$line_info[$j]})) {
	  push @sample_keep_inds, 0;
	  next;
	}

	push @sample_keep_inds, 1;
	if ($is_samp_snptest) {
	  print SAMP "$line_info[$j] $line_info[$j] 0\n";
	}
	else {
	  print SAMP "$line_info[$j]\n";
	}
      }
      close(SAMP);

      last; # done reading header
    }
  }

  # process SNP information
  my $a = 0;
  my $n_filtered_snps = 0;
  my $n_multiallelic_snps = 0;
  my $n_monomorphic_snps = 0;
  my $n_non_snp_variants = 0;
  my $n_non_indel_variants = 0;
  my $n_non_sv_variants = 0;
  my $is_chrom = 0;  # have we reached the chromosome of interest yet?
  my $gt_like_idx = -1;
  my (@format_info, @gtype_info, @gt_like_info, @alleles);
  while (defined ($line = <VCF>)) {
    chomp $line;
    @line_info = split /\t/, $line;
    #print STDERR "|$line_info[0]| -- ".scalar(@line_info)."\n";

    if (++$a % 1000 == 0) {
      print STDERR "$line_info[0] -- $line_info[$pos_tag_idx]\n";
    }

    $is_chrom = 1 if ($line_info[0] eq $chrom);
    if ($chrom !~ m/all$/i && $line_info[0] ne $chrom) {
      last if ($is_chrom);  # assumes that all VCF rows for this chrom are contiguous
      next;
    }

    next if ($line_info[$pos_tag_idx] < $start_pos);
    last if ($line_info[$pos_tag_idx] > $end_pos);

    if ($line_info[$filter_tag_idx] ne "PASS" && $line_info[$filter_tag_idx] ne ".") {
      ++$n_filtered_snps;
      next;
    }

    if ($line_info[$alt_tag_idx] =~ ',') {
      ++$n_multiallelic_snps;
      next;
    }

    if ($line_info[$alt_tag_idx] eq '.' && $no_mono) {
      ++$n_monomorphic_snps;
      next;
    }

    # SNPs are defined as variants that have two (and only two) single-base alleles
    if ($snps_only &&
	!(exists($snp_alleles{$line_info[$ref_tag_idx]}) && exists($snp_alleles{$line_info[$alt_tag_idx]})) ) {
      ++$n_non_snp_variants;
      next;
    }

    # INDELs are defined as variants that include one single-base allele and one allele with multiple bases;
    # this follows the convention used in VCF files from the 1000 Genomes Project
    if ($indels_only &&
	!((exists($snp_alleles{$line_info[$ref_tag_idx]}) && length($line_info[$alt_tag_idx]) > 1) ||
	  (length($line_info[$ref_tag_idx]) > 1 && exists($snp_alleles{$line_info[$alt_tag_idx]})) )) {
      ++$n_non_indel_variants;
      next;
    }

    # SVs are defined as variants for which neither allele is a single-base allele \in {A,C,G,T}, *unless*
    # one of the alleles is '-'; this follows the convention used in VCF files from the 1000 Genomes Project
    if ($svs_only &&
	!($line_info[$ref_tag_idx] eq '-' || $line_info[$alt_tag_idx] eq '-' ||
	  (length($line_info[$ref_tag_idx]) > 1 && length($line_info[$alt_tag_idx]) > 1)) ) {
      ++$n_non_sv_variants;
      next;
    }

    # change coding for deletion alleles
    my $ref_allele = $line_info[$ref_tag_idx];
    my $alt_allele = $line_info[$alt_tag_idx];
    $ref_allele = "-" if ($ref_allele eq "<DEL>");
    $alt_allele = "-" if ($alt_allele eq "<DEL>");

    # this SNP passed filtering; print header info and process the genotypes
    print GEN "$line_info[0] $line_info[$id_tag_idx] $line_info[$pos_tag_idx] $ref_allele $alt_allele";

    # if we want to parse and print genotype likelihoods rather than hard calls,
    # use the format string to determine which field contains this information
    $gt_like_idx = -1;
    if ($is_gt_like) {
      @format_info = split /\:/, $line_info[$format_tag_idx];
      for my $i (0..$#format_info) {
	if ($format_info[$i] eq "GL") {
	  $gt_like_idx = $i;
	  last;
	}
      }
    }

    my $j = -1;
    for my $i (($format_tag_idx+1)..$#line_info) {
      ++$j;
      next if (!$sample_keep_inds[$j]);

      @gtype_info = split /\:/, $line_info[$i];
      if ($is_gt_like && $gt_like_idx > -1) { # print genotype likelihoods
	@gt_like_info = split /,/, $gtype_info[$gt_like_idx];
	die $! if (scalar(@gt_like_info) != 3);

	# back-convert from log10 to probability scale
	my $tot = 0;
	for my $k (0..2) {
	  $gt_like_info[$k] = 10**$gt_like_info[$k];
	  $tot += $gt_like_info[$k];
	}
	$tot = 1.0 if ($tot == 0);

	# normalize probs and print
	for my $k (0..2) {
	  $gt_like_info[$k] /= $tot;
	  print GEN " ".sprintf("%.4f",$gt_like_info[$k]);
	}
      }
      else { # print hard-call genotypes
	if ($is_special_gt && exists($special_gt_to_probs{$gtype_info[0]})) {
	  print GEN " $special_gt_to_probs{$gtype_info[0]}";
	}
	elsif (exists($unphased_gt_to_probs{$gtype_info[0]})) {
	  print GEN " $unphased_gt_to_probs{$gtype_info[0]}";
	}
	elsif (exists($phased_gt_to_probs{$gtype_info[0]})) {
	  print GEN " $phased_gt_to_probs{$gtype_info[0]}";
	}
	else {
	  die "\nERROR: Genotype |$gtype_info[0]| (i=$i) at site $line_info[$pos_tag_idx] is not in accepted format.\n\n";
	}
      }
    }
    print GEN "\n";
  }

  # print number of SNPs removed by FILTER tag
  print STDERR "\n--$n_filtered_snps SNPs were removed by the FILTER tag.\n";
  print STDERR "\n--$n_multiallelic_snps SNPs were removed for having more than one allele in the ALT column.\n";
  print STDERR "\n--$n_monomorphic_snps SNPs were removed for having no defined allele in the ALT column.\n" if ($no_mono);
  print STDERR "\n--$n_non_snp_variants variants were removed for having non-SNP alleles.\n" if ($snps_only);
  print STDERR "\n--$n_non_indel_variants variants were removed for having non-INDEL alleles.\n" if ($indels_only);
  print STDERR "\n--$n_non_sv_variants variants were removed for having non-SV alleles.\n" if ($svs_only);
  #print STDERR "\n--$n_low_maf_snps variants were removed for having MAFs below the -min_maf threshold.\n" if ($n_low_maf_snps > 0);
  print STDERR "\n";

  close(VCF);
  close(GEN);
}


# parse the command line
sub ParseCommandLine()
{
  my $usage = "\nvcf2impute_gen.pl\n";
  $usage .= "  [-vcf]           VCF file to process; can be phased or unphased\n";
  $usage .= "  [-gen]           name of new file to print; automatically gzipped\n";
  $usage .= "  <-chr>           chromosome to include in output files, in (chr)[1-22,X]\n";
  $usage .= "  <-start>         first position to include in output files\n";
  $usage .= "  <-end>           last position to include in output files\n";
  $usage .= "  <-samp_include>  file containing a list of samples to include in gen file\n";
  $usage .= "  <-samp_exclude>  file containing a list of samples to exclude from gen file\n";
  $usage .= "  <-samp_snptest>  flag: print sample file in SNPTEST format\n";
  $usage .= "  <-special_gt>    flag: allow 'special' genotype coding in VCF\n";
  $usage .= "  <-gt_like>       flag: print genotype likelihoods rather than hard calls\n";
  $usage .= "  <-no_mono>       flag: omit monomorphic sites from processed files\n";
  $usage .= "  <-snps_only>     flag: omit any sites that are not biallelic SNPs\n";
  $usage .= "  <-indels_only>   flag: omit any sites that are not biallelic INDELs\n";
  $usage .= "  <-svs_only>      flag: omit any sites that are not biallelic SVs\n";
  $usage .= "\n";
  $usage .= "  (args in square brackets required; args in pointy brackets optional)\n";
  $usage .= "\n";

  # default values
  $arg{-vcf} = "";
  $arg{-gen} = "";
  $arg{-chr} = "All";
  $arg{-start} = 0;
  $arg{-end} = 1e9;
  $arg{-samp_include} = "";
  $arg{-samp_exclude} = "";
  $arg{-samp_snptest} = 0;  # boolean variable; false unless flag provided
  $arg{-special_gt} = 0;  # boolean variable; false unless flag provided
  $arg{-gt_like} = 0;  # boolean variable; false unless flag provided
  $arg{-no_mono} = 0;  # boolean variable; false unless flag provided
  $arg{-snps_only} = 0;  # boolean variable; false unless flag provided
  $arg{-indels_only} = 0;  # boolean variable; false unless flag provided
  $arg{-svs_only} = 0;  # boolean variable; false unless flag provided

  # parse the command line
  for (my $i = 0; $i <= $#ARGV; ++$i) {
    if ($ARGV[$i] =~ /^-/) {
      if ($ARGV[$i] eq "-samp_snptest" || $ARGV[$i] eq "-special_gt" || $ARGV[$i] eq "-gt_like" || $ARGV[$i] eq "-no_mono" ||
	  $ARGV[$i] eq "-snps_only" || $ARGV[$i] eq "-indels_only" || $ARGV[$i] eq "-svs_only") {
	$arg{$ARGV[$i]} = 1;  # set boolean to 'true'
      }
      else {
	$arg{$ARGV[$i]} = $ARGV[$i+1];
      }
    }
  }

  die ($usage) if ($arg{-vcf} eq "");
  die ($usage) if ($arg{-gen} eq "");
}
