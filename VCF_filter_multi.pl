#!/usr/bin/env perl
#Author: Yulong Li <liyulong12@mails.ucas.ac.cn>
use strict;
use warnings;
use threads;
use threads::shared;
use List::MoreUtils qw{uniq};
use Getopt::Long qw(:config no_ignore_case bundling);
use File::Spec;
use File::Basename;
my $abs_path = File::Spec->rel2abs( __FILE__ );
my @args = @ARGV;
my @out :shared = ();
my ($cnt, $snp_total) :shared;
my (@header, $pops, @order, $in_fh, $out_fh, %samples, $tot_indiv, $num_pops);
my ($cmd, $sort, $infile, $outfile, $popmap, $minDP, $maxDP, $het, $fis, $hwe);
my ($cov, $tot_cov, $minGQ, $minQ, $l_maf, $g_maf, $global, $num_threads, $help);
my ($totDP,$avgDP, $biallelic, $polymorphic, $rmIndel, $poly_thred);
my $vcf_pre_header = '';
# Initialize
$num_threads = 1;
$cnt         = 0;
$snp_total   = 0;
$biallelic   = 1;
$polymorphic = 1;
$rmIndel     = 1;

my $x ='#';
GetOptions (
    "in=s"          => \$infile,
    "out=s"         => \$outfile,
    "Popmap=s"      => \$popmap,
    "minDP=s"       => \$minDP,
    "MaxDP=s"       => \$maxDP,
    "Het=f"         => \$het,
    "Fis=f"         => \$fis,
    "phwe=f"        => \$hwe,
    "GQ=i"          => \$minGQ,
    "Q=i"           => \$minQ,
    "localMAF=f"    => \$l_maf,
    "globalMAF=f"   => \$g_maf,
    "filter:1"      => \$global,
    "sort:1"        => \$sort,
    "Threads=i"     => \$num_threads,
    "coverage=f"    => \$cov,
    "Cov=f"         => \$tot_cov,
    "totDP=s"       => \$totDP,
    "avgDP=s"       => \$avgDP,
    "reNb:0"        => \$biallelic,
    "reMo:0"        => \$polymorphic,
    "reIn:0"        => \$rmIndel,
	"rePo=f"        => \$poly_thred,
    "help:1"        => \$help
) or die ("Error in command line arguments\n");

my $usage = 
'   
    Options [* required, case sensitive]:
    
    --i *input vcf file.
    --o *output file.
    --P *PopMap file.
    --c  coverage for each population.
    --C  total genotyping rate.
    --m  min depth of each individual.
    --M  max depth of each individual.
    --H  max Ho of each population.
    --F  abs Fis value.
    --G  min Genotype quality.
    --Q  min quality for each SNP.
    --g  Global MAF.
    --l  Local MAF.
    --p  p-value for hwe.
    --f  apply the filter (Ho and Fis) on each site instead of each population.
    --T  number of threads.
    --s  sort the final vcf file.
    --h  help.
    --totDP total Depth range for each site, [min:max].
    --avgDP average Depth range for each site, [min:max].
    --reNb  retain non-bi-allelic sites.
    --reMo  retain monomorphic sites.
    --reIn  retain indels.
    --rePo  threshold for polymorphic loci of each populations.
	
';
    
die "$usage\n" if $help or !($infile && $outfile && $popmap);
my @file_stat = stat ($infile);
my $max_per_batch = 999;
if (substr($infile, -2) eq 'gz') {
    open($in_fh, "gunzip -c $infile|") or die "$!";
    $max_per_batch = $file_stat[7]/1e6 > 100 ? 9999 : 999; # 100M
    } 
elsif (substr($infile, -2) eq '7z') {
    open($in_fh, "7zcat.sh $infile|") or die "$!";
    $max_per_batch = $file_stat[7]/1e6 > 100 ? 9999 : 999; # 100M
    }
elsif (substr($infile, -3) eq 'bcf') {
    open($in_fh, "bcftools convert --threads 8 -Ov $infile |") or die "$!";
    $max_per_batch = $file_stat[7]/1e6 > 100 ? 9999 : 999; # 100M
    }
else {
    open($in_fh, "$infile") or die "$!";
    $max_per_batch = $file_stat[7]/1e6 > 1000 ? 9999 : 999; # 1G
}

###### Parameters ######
no warnings;
my $H = $global? "Global Ho" : "MaxHo";
my $F = $global? "Global Abs Fis" : "Abs Fis";
my $para  = "CMD: $abs_path ".join(' ', @args);
   $para .= "\n\nParameters used:\n";
   $para .= sprintf "\t%-20s : %-10s\n", 'File name', basename($infile);
   $para .= sprintf "\t%-20s : %-10s\n", 'Out name', basename($outfile);
   $para .= sprintf "\t%-20s : %-10s\n", 'Local coverage', $cov if $cov;
   $para .= sprintf "\t%-20s : %-10s\n", 'Total coverage', $tot_cov if $tot_cov;
   $para .= sprintf "\t%-20s : %-10s\n", 'MinGQ', $minGQ if $minGQ;
   $para .= sprintf "\t%-20s : %-10s\n", 'MinQ', $minQ if $minQ;
   $para .= sprintf "\t%-20s : %-10s\n", 'minDepth', $minDP if $minDP;
   $para .= sprintf "\t%-20s : %-10s\n", 'MaxDP', $maxDP if $maxDP;
   $para .= sprintf "\t%-20s : %-10s\n", $H, $het if $het;
   $para .= sprintf "\t%-20s : %-10s\n", $F, $fis if $fis;
   $para .= sprintf "\t%-20s : %-10s\n", 'p-value for hwe test', $hwe if $hwe;
   $para .= sprintf "\t%-20s : %-10s\n", 'Total Depth', $totDP if $totDP;
   $para .= sprintf "\t%-20s : %-10s\n", 'Average Depth', $avgDP if $avgDP;
   $para .= sprintf "\t%-20s : %-10s\n", 'Global MAF', $g_maf if $g_maf;
   $para .= sprintf "\t%-20s : %-10s\n", 'Local MAF', $l_maf if $l_maf;
   $para .= sprintf "\t%-20s : %-10s\n", 'Local polymorphic', $poly_thred if defined $poly_thred;
   $para .= sprintf "\t%-20s\n", 'Only keep polymorphic sites' if $polymorphic;
   $para .= sprintf "\t%-20s\n", 'Only keep bi-allelic sites' if $biallelic;
   $para .= sprintf "\t%-20s\n", 'Remove Indels' if $rmIndel;
use warnings;
open($out_fh, ">$outfile.log") or die "$!";
print $out_fh $para;
close $out_fh;
print STDERR "\n=================================================================================\n";
print STDERR "Using $num_threads threads, $max_per_batch per batch...\n\n";
print STDERR "$para\n";

if (substr($outfile, -2) eq 'gz') {
    # use bgzip.
    open($out_fh, "|bgzip -c >$outfile") or die "$!";
    } else {
    open($out_fh, ">$outfile") or die "$!";
}

my $head = 0;
while (<$in_fh>) {
    #Just output number of SNPs if nothing need to be filtered.
    unless ($cov || $minGQ || $minDP || $maxDP || $het || $fis
    || $g_maf || $l_maf || $minQ || $tot_cov || $hwe || $polymorphic || $biallelic || $rmIndel) {
    
        print $out_fh $_;
        next if /^##|^$/;
        if (/^#/) {
        chomp;
        my @tmp = split(/\t/);
        $tot_indiv = @tmp - 9;
        next;
        }
        $cnt++;
        $snp_total = $cnt;
        next;
    }
    next if /^$/;
    if (/^##/) {
        $vcf_pre_header .= $_;
        next;
    }
    next if $head > 0 && /^#/; 
    if ($head ==0 && /^#/) {
        chomp;
        push @header, split(/\t/);
        for (my $i=9; $i<=$#header; $i++) {
            $samples{$header[$i]} = $i;
        }
        
        ##################
        
        #Read PopMap file;
        
        #[Popmap]\t[file]:
        #indiv1\tpop1
        #indiv2\tpop2
        #indiv3\tpop1 ...
        
        ##################
        
        open(my $pop_fh, "$popmap") or die "No PopMap file!";
        while (<$pop_fh>) {
            next if /^#|^$/;
            $_ =~ s/[\r\n]|^\s+|\s+$//g;
            my @part = split;
            my $cpop = $part[1];
            my $cind = $part[0];
            die "No such individual $cind in VCF\n" if not defined $samples{$cind};
            push @order, $cpop;
            push @{$pops->{$cpop}}, $samples{$cind}; #Pop name => @indiv rank. 
            
        }
        close $pop_fh;
        #die "Header is wrong!" if @order != (@header -9);
        $tot_indiv = @order; 
        @order     = uniq @order;
        $num_pops  = scalar(@order);
        ###### Print header ######
        print STDERR "PopMap($num_pops):\n";
        print $out_fh $vcf_pre_header;
        print $out_fh join("\t", @header[0..8]);
        foreach my $pop (@order) {
            my @indivs     = @header[@{$pops->{$pop}}];
            my $num_indivs = scalar(@indivs);
            print $out_fh "\t", join("\t", @indivs);
            printf STDERR "%-s(%-d): %-s;\n", $pop, $num_indivs, join(", ", @indivs);
        }
        print STDERR "\n";
        print $out_fh "\n";
        $head++;
        next;
        
    }
    
    my @batch = ($_);
    my @self  = ();
    for (1..$num_threads) {
        for (1..$max_per_batch) {
            last if eof;
            my $vcf = <$in_fh>;
            push @batch, $vcf; 
        }
        my $t = threads->create(\&filter, \@batch);
        push @self, $t;
        undef @batch;
        last if eof;
    }
    map { $_->join();} @self;
    print STDERR "Writing file...$cnt of $snp_total\r" if $cnt;
    print $out_fh @out;
    undef @out;
}
close $out_fh;
print STDERR "\ndone.\n";

if ($sort) {
    print STDERR "Sorting file...";
    if (substr($outfile, -2) eq 'gz') {
        # use bgzip.
        `(zgrep '^#' $outfile; zgrep -v '^#' $outfile | sort -T ./ --parallel=$num_threads -k1,1V -k2,2n ) | bgzip --threads 4 -c >$outfile.1 && mv $outfile.1 $outfile`;
    } else {
        `(grep '^#' $outfile; grep -v '^#' $outfile | sort -T ./ --parallel=$num_threads -k1,1V -k2,2n) >$outfile.1 && mv $outfile.1 $outfile`;
    }
    print STDERR "done.\n";
}
open($out_fh, ">>$outfile.log") or die "$!";
print $out_fh "\nFinal retained\t: $cnt SNPs of $tot_indiv individuals\nTotal $snp_total SNPs\n";
close $out_fh;
print STDERR "=================================================================================";
print STDERR "\nFinal retained\t: $cnt SNPs of $tot_indiv individuals\nTotal $snp_total SNPs\n";

sub read_file {
    my $fh    = shift;
    my $line  = shift;
    my @lines = ();
    map {push @lines, <$fh>} 1..$line;
    return \@lines;
}

sub calc_stat {

    my ($gts, $ref, $alt) = @_;
    my $allele->{$ref} = 0;
       $allele->{$alt} = 0;
       $allele->{'het'}= 0;
    my $obs_hets = 0;
    my $obs_hom1 = 0;
    my $obs_hom2 = 0;
    
    foreach my $gt (@$gts) {
    
        # unphased or phased.
        if ($gt eq '0/0' || $gt eq '0|0') {$allele->{$ref} += 2; $obs_hom1++;}
        if ($gt eq '1/1' || $gt eq '1|1') {$allele->{$alt} += 2; $obs_hom2++;}
        if ($gt eq '0/1' || $gt eq '0|1' || $gt eq '1|0' || $gt eq '1/0') {
            $allele->{$ref}++; 
            $allele->{$alt}++;
            $allele->{'het'}++;
            $obs_hets++;
        }
    }
    my $n    = ($allele->{$ref} + $allele->{$alt}) / 2;
    my $p    = $allele->{$ref} / (2 * $n);
    my $Ho   = $allele->{'het'} / $n;
    my $He   = 2 * $p * (1 - $p);
    my $Fis  = $He > 0 ? (1 - $Ho/$He) : 'NAN';
    my $flag = $p < 0.5 ? $ref : $alt;
    my $maf  = $p < 0.5 ? $p : (1 - $p);
    my $p_hwe= 1; $p_hwe = hwe($obs_hets, $obs_hom1, $obs_hom2) if defined($hwe);
    return($Ho, $He, $Fis, $maf, $flag, $p_hwe);
    
}

sub hwe {
# snphwe.pl: A Perl implementation of the fast exact Hardy-Weinberg Equilibrium 
# test for SNPs as described in Wigginton, et al. (2005). 
#
# Copyright 2010 Joshua Randall
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Author:
#    This software was written Joshua C. Randall <jcrandall@alum.mit.edu>
#    and is a port to Perl of algorithms implemented by others in C and R.
#
# Attribution:
#    This software is based entirely on the C and R implementations of the 
#    algorithms for exact HWE tests as described in Wigginton, et al. (2005), 
#    which were originally written by Jan Wigginton and released into the 
#    public domain.  C, R, and Fortran implementations of these algorithms 
#    are available for download at:
#       http://www.sph.umich.edu/csg/abecasis/Exact/
#
# Citation: 
#    This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
#    Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of 
#    Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76(5): 887 - 893.
#    Please cite this work when using this code.
#
# Usage: 
#    This software is a Perl library, intended to be used within other programs.
#    To use this library directly from the command-line, you can run Perl with a 
#    one-line program such as:
#
#       perl -e 'require "snphwe.pl"; print(snphwe(@ARGV))' 57 14 50
#
#    Where the three numbers at the end are the observed counts of the three 
#    genotypes: first the heterozygote count, then one of the homozygote genotype 
#    counts, and finally the other homozygote genotype count, in that order.  
#
#    The example above, which would be for 57 Aa, 14 aa, and 50 AA, should print 
#    the resulting P-value, which in this case is 0.842279756570793, to the 
#    standard output.
#
# Note:
#    Code for the alternate P-value calculation based on p_hi/p_lo that was 
#    included in the Wigginton, et al. C and R implementations (but was 
#    disabled) has been included here, but has not been tested.  It is 
#    therefore commented out.  If you wish to make use of this code, please 
#    verify it functions as desired.
#
    my $obs_hets = shift;
    my $obs_hom1 = shift;
    my $obs_hom2 = shift;

    if($obs_hom1 < 0 || $obs_hom2 < 0 || $obs_hets <0) {
	return(-1);
    }

    # rare homozygotes
    my $obs_homr;

    # common homozygotes
    my $obs_homc;
    if($obs_hom1 < $obs_hom2) {
	$obs_homr = $obs_hom1;
	$obs_homc = $obs_hom2;
    } else {
	$obs_homr = $obs_hom2;
	$obs_homc = $obs_hom1;
    }

    # number of rare allele copies
    my $rare_copies = 2 * $obs_homr + $obs_hets;

    # total number of genotypes
    my $genotypes = $obs_homr + $obs_homc + $obs_hets;

    if($genotypes <= 0) {
	return(-1);
    }
    
    # Initialize probability array
    my @het_probs;
    for(my $i=0; $i<=$rare_copies; $i++) {
	$het_probs[$i] = 0.0;
    }

    # start at midpoint
    my $mid = int($rare_copies * (2 * $genotypes - $rare_copies) / (2 * $genotypes));

    # check to ensure that midpoint and rare alleles have same parity
    if(($rare_copies & 1) ^ ($mid & 1)) {
	$mid++;
    }
    
    my $curr_hets = $mid;
    my $curr_homr = ($rare_copies - $mid) / 2;
    my $curr_homc = $genotypes - $curr_hets - $curr_homr;

    $het_probs[$mid] = 1.0;
    my $sum = $het_probs[$mid];
    for($curr_hets = $mid; $curr_hets > 1; $curr_hets -= 2) {
	$het_probs[$curr_hets - 2] = $het_probs[$curr_hets] * $curr_hets * ($curr_hets - 1.0) / (4.0 * ($curr_homr + 1.0) * ($curr_homc + 1.0));
	$sum += $het_probs[$curr_hets - 2];

	# 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
	$curr_homr++;
	$curr_homc++;
    }

    $curr_hets = $mid;
    $curr_homr = ($rare_copies - $mid) / 2;
    $curr_homc = $genotypes - $curr_hets - $curr_homr;
    for($curr_hets = $mid; $curr_hets <= $rare_copies - 2; $curr_hets += 2) {
	$het_probs[$curr_hets + 2] = $het_probs[$curr_hets] * 4.0 * $curr_homr * $curr_homc / (($curr_hets + 2.0) * ($curr_hets + 1.0));
	$sum += $het_probs[$curr_hets + 2];
	
	# add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
	$curr_homr--;
	$curr_homc--;
    }

    for(my $i=0; $i<=$rare_copies; $i++) {
	$het_probs[$i] /= $sum;
    }

    # alternate p-value calculation for p_hi/p_lo
#    my $p_hi = $het_probs[$obs_hets];
#    for(my $i=$obs_hets+1; $i<=$rare_copies; $i++) {
#	$p_hi += $het_probs[$i];
#    }
#    
#    my $p_lo = $het_probs[$obs_hets];
#    for(my $i=$obs_hets-1; $i>=0; $i--) {
#	$p_lo += $het_probs[$i];
#    }
#
#    my $p_hi_lo;
#    if($p_hi < $p_lo) {
#	$p_hi_lo = 2 * $p_hi;
#    } else {
#	$p_hi_lo = 2 * $p_lo;
#    }

    # Initialise P-value 
    my $p_hwe = 0.0;

    # P-value calculation for p_hwe
    for(my $i = 0; $i <= $rare_copies; $i++) {
	if($het_probs[$i] > $het_probs[$obs_hets]) {
	    next;
	}
	$p_hwe += $het_probs[$i];
    }
    
    if($p_hwe > 1) {
	$p_hwe = 1.0;
    }

    return($p_hwe);
}

sub filter {

    my $vcf = shift;
    foreach my $vcfline (@$vcf) {
    #################################################################
    
    #VCF format:
    #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Indiv_genotypes ...
    #Genotypes:GT:PL:DP:SP:GQ ...
    
    #Filter order: biallelic -> Q -> depth -> GQ -> coverage -> local Ho or Fis or MAF -> hwe
    #-> missing rate -> global Ho or Fis or MAF.
    
    #################################################################
    next if $vcfline =~ /^#/; 
    chomp $vcfline;
    
    if (1) {
        lock($snp_total);
        $snp_total++; # Number of total SNPs.
    }
    
    my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @genos) = split(/\t/, $vcfline);
    my $fmt;
    my @filtered         = ();
    my @gts_total        = ();
    my $tot_miss         = 0; # Missing rate for a site.
    my $cnt_Fis          = 0;
    my $cnt_Ho           = 0;
	my $cnt_poly         = 0; # not a polymorphic loci for one population.
    my $cnt_hwe          = 0;
    my $cnt_gmaf         = 0;
    my $cnt_lmaf->{$ref} = 0;
       $cnt_lmaf->{$alt} = 0;
    my $mono             = 0;
    my $tDP              = 0;
    next if not ($filter eq "." || $filter eq "PASS");
    next if ($biallelic && $alt =~ /\,|\./);        # Skip non-biallelic loci.
    next if ($rmIndel && ($ref =~ /\w{2,}/ || $alt =~ /\w{2,}/)); # remove indel.
    if(defined($minQ)) {
        # Minmum quality.
        if($qual eq '.') {
            print STDERR "No Quality value in VCF lines.\n";
        } else {
            next if ($qual < $minQ);
        }
    }
    my $recode = 0;                                 # Populations number of non-enough coverage. 
    my @formats = split(/:/, $format);              # Order of format:
    map { $fmt->{$formats[$_]} = $_;} 0..$#formats; # Geno => Order.
    
    #Iteration each population.
    foreach my $name (@order) {

        my $total     = @{$pops->{$name}};
        my $miss      = 0;
        my $cov_ratio = 1;
        my @gts       = (); #Genotypes for each population.
        
        foreach my $rank(@{$pops->{$name}}) {
            #each individual of every population.
            my $geno0= $genos[$rank-9]; # Genotype of original.
            my @geno = split(/:/, $geno0);
            if (not defined $geno[$fmt->{'GT'}]) { die "GT field is requred!\n";}
            #if (not defined $geno[$fmt->{'DP'}]) { die "DP field is requred!\n";}
            my $GT   = $geno[$fmt->{'GT'}];
            
            #
            # Missing is skipped.
            #
            if ($GT eq './.' || $GT eq '.' || $GT eq '.|.') {
            
                push (@filtered, $geno0);
                $miss++;
                $tot_miss++;
                next;
            }
            
            my $DP   = $geno[$fmt->{'DP'}] if $fmt->{'DP'};
            my $GQ   = $geno[$fmt->{'GQ'}] if $fmt->{'GQ'};
            $tDP += $DP; # total depth.
            
            #Other filter...
            ###### Depth ######
            if ((defined($minDP) && $DP < $minDP) or (defined($maxDP) && $DP > $maxDP)) {
                
                $geno[$fmt->{'GT'}] = './.';
                $miss++;
                $tot_miss++;
                push @filtered, join(":", @geno);
                next;
            }
            
            ###### GQ ######
            if (defined($minGQ) && $GQ < $minGQ) {
                $geno[$fmt->{'GT'}] = './.';
                $miss++;
                $tot_miss++;
                push @filtered, join(":", @geno);
                next;
            }
            
            push (@filtered, $geno0);
            push (@gts, $GT); #For calculate statistics, Skip ./.
            push (@gts_total, $GT); #Global genotypes, Skip ./.
            
        }
        
        ###### Coverage ######
        my $l_N = $total - $miss; # local num of individuals.
        if (defined $cov) {
            
            $cov_ratio = $cov > 1 ? $l_N : $l_N/$total;
            $recode++ if $cov_ratio < $cov;
            last if $cov_ratio < $cov; # next site.
        } else {
            next if $l_N == 0; # next population.
        }
        ###### Local ######
        if ((!$global && (defined($het) or defined ($fis))) or defined($l_maf) or defined($hwe)) {
        
            my ($Ho, $He, $Fis, $maf, $flag, $p_hwe) = calc_stat(\@gts, $ref, $alt); #Each population.
            my $mac = sprintf "%.1f", 2*$l_N*$maf; # try to not lose precision.
			$maf = ($l_maf && $l_maf > 1) ? $mac : $maf; # allele count or freq.
            if (!$global) {
                ###### Het ######
                if (defined($het) && $Ho > $het) {
                    $cnt_Ho++;
                    last;
                }
            
                ###### Fis ######
                if (defined($fis) && abs($Fis) > $fis) {
                    $cnt_Fis++;
                    last;
                }
                ###### HWE ######
                $cnt_hwe++ if (defined ($hwe) && $p_hwe < $hwe);
            }
            ###### MAF ######
            $cnt_lmaf->{$flag}++ if $l_maf && $maf < $l_maf;

        }
		# polymorphic loci of each population.
		if (defined $poly_thred) {
			my ($Ho, $He, $Fis, $maf, $flag, $p_hwe) = calc_stat(\@gts, $ref, $alt);
			my $mac = sprintf "%.1f", 2*$l_N*$maf;
			my $maf_= $poly_thred >= 1 ? $mac : $maf;
			if ($maf_ < $poly_thred) {
				$cnt_poly++;
				last;
			}
		}
        
    }
    
    next if ($recode  > 0 || $cnt_Ho  > 0 || $cnt_Fis > 0 || $cnt_poly > 0);
    next if (defined($hwe) && $cnt_hwe >= $num_pops);
    
    ###### Global missing rate #######
    my $g_N = $tot_indiv - $tot_miss; # global num of individuals.
    if (defined $tot_cov) {
        my $cov_ratio = $tot_cov > 1 ? $g_N : $g_N/$tot_indiv;
        next if $cov_ratio < $tot_cov;
    }
    ###### total depth and average depth ######
    if (defined $totDP) {
        if ($totDP =~ /(\d+):(\d+)/) {
            my $low  = $1;
            my $high = $2;
            next if ($tDP > $high || $tDP < $low);
        } else {
            next if $tDP < $totDP;
        }
    }
    
    if (defined $avgDP) {
        my $aDP = $tDP / $g_N; # average depth.
        if ($avgDP =~ /(\d+):(\d+)/) {
            my $low  = $1;
            my $high = $2;
            next if ( $aDP > $high || $aDP < $low);
        } else {
            next if $aDP < $totDP;
        }
    
    }
    
    ###### Global ######
    if (defined($global) or defined($g_maf)) {
        
        my ($Ho, $He, $Fis, $maf, $flag, $p_hwe) = calc_stat(\@gts_total, $ref, $alt);
        $maf = ($g_maf && $g_maf > 1) ? 2*$g_N*$maf : $maf;
        if ($global) {
            ###### Het ######
            if (defined($het) && $Ho > $het) {
                $cnt_Ho++;
                next;
            }
            ###### Fis ######
            if (defined($fis) && abs($Fis) > $fis) {
                $cnt_Fis++;
                next;
            }
            ###### HWE ######
            if (defined($hwe) && $p_hwe < $hwe) {$cnt_hwe++;next;}
        }
        
        ###### MAF ######
        $cnt_gmaf++ if defined($g_maf) && $maf < $g_maf;
        if (not defined $l_maf) {
            $cnt_lmaf->{$ref} = $num_pops;
            $cnt_lmaf->{$alt} = $num_pops;
        }
        
    }
       
    next if ($cnt_Ho   > 0 || $cnt_Fis  > 0);
    next if $cnt_gmaf > 0 && ($cnt_lmaf->{$ref} == $num_pops or $cnt_lmaf->{$alt} == $num_pops);
    ## check if is snp.
    next if ($polymorphic && scalar(uniq @gts_total) == 1);
    
    ###### Print each snp sites ######
    lock(@out);
    lock($cnt);
    push @out, join("\t", $chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @filtered) . "\n";
    $cnt++; # SNP number.
    #printf STDERR "Filtering: %-1s [Retained loci: %d]\r", $x, $cnt;
    #printf STDERR "Filtering: retained %d of total %d loci, %d individuals of %d populations.\r", $cnt, $snp_total, $tot_indiv, $num_pops;
    }
}