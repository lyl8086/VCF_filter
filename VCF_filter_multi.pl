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
my ($cmd, $sort, $infile, $outfile, $popmap, $minDP, $maxDP, $het, $fis);
my ($cov, $tot_cov, $minGQ, $minQ, $l_maf, $g_maf, $global, $num_threads, $help);
my $x ='#';
GetOptions (
    "in=s"          => \$infile,
    "out=s"         => \$outfile,
    "Popmap=s"      => \$popmap,
    "minDP=s"       => \$minDP,
    "MaxDP=s"       => \$maxDP,
    "Het=f"         => \$het,
    "Fis=f"         => \$fis,
    "GQ=i"          => \$minGQ,
    "Q=i"           => \$minQ,
    "localMAF=f"    => \$l_maf,
    "globalMAF=f"   => \$g_maf,
    "filter:1"      => \$global,
    "sort:1"        => \$sort,
    "threads=i"     => \$num_threads,
    "coverage=f"    => \$cov,
    "Cov=f"         => \$tot_cov,
    "help:1"        => \$help
) or die ("Error in command line arguments\n");

my $usage = 
'   
    Options [* required]:
    
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
    --f  apply the filter (Ho and Fis) on each site instead of each population.
    --t  number of threads.
    --s  sort the final vcf file.
    --h  help.
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
   $para .= sprintf "\t%-15s : %-10s\n", 'File name', basename($infile);
   $para .= sprintf "\t%-15s : %-10s\n", 'Out name', basename($outfile);
   $para .= sprintf "\t%-15s : %-10s\n", 'Local coverage', $cov if $cov;
   $para .= sprintf "\t%-15s : %-10s\n", 'Total coverage', $tot_cov if $tot_cov;
   $para .= sprintf "\t%-15s : %-10s\n", 'MinGQ', $minGQ if $minGQ;
   $para .= sprintf "\t%-15s : %-10s\n", 'MinQ', $minQ if $minQ;
   $para .= sprintf "\t%-15s : %-10s\n", 'minDepth', $minDP if $minDP;
   $para .= sprintf "\t%-15s : %-10s\n", 'MaxDP', $maxDP if $maxDP;
   $para .= sprintf "\t%-15s : %-10s\n", $H, $het if $het;
   $para .= sprintf "\t%-15s : %-10s\n", $F, $fis if $fis;
   $para .= sprintf "\t%-15s : %-10s\n", 'Global MAF', $g_maf if $g_maf;
   $para .= sprintf "\t%-15s : %-10s\n", 'Local MAF', $l_maf if $l_maf;
use warnings;
open($out_fh, ">$outfile.log") or die "$!";
print $out_fh $para;
close $out_fh;
print STDERR "\n=================================================================================\n";
print STDERR "Using $num_threads threads, $max_per_batch per batch...\n\n";
print STDERR "$para\n";

if (substr($outfile, -2) eq 'gz') {
    open($out_fh, "|gzip >$outfile") or die "$!";
    } else {
    open($out_fh, ">$outfile") or die "$!";
}

my $head = 0;
while (<$in_fh>) {
    #Just output number of SNPs if nothing need to be filtered.
    unless ($cov || $minGQ || $minDP || $maxDP || $het || $fis
    || $g_maf || $l_maf || $minQ || $tot_cov) {
    
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

    next if /^##|^$/;
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
        
        open(my $pop, "$popmap") or die "No PopMap file!";
        while (<$pop>) {
            next if /^#|^$/;
            $_ =~ s/[\r\n]|^\s+|\s+$//g;
            my @part = split;
            push @order, $part[1];
            push @{$pops->{$part[1]}}, $samples{$part[0]}; #Pop name => @indiv rank. 
            
        }
        close $pop;
        die "Header is wrong!" if @order != (@header -9);
        $tot_indiv = @order; 
        @order     = uniq @order;
        $num_pops  = scalar(@order);
        ###### Print header ######
        print STDERR "PopMap($num_pops):\n";
        print $out_fh join("\t", @header[0..8]);
        foreach $pop (@order) {
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
        `(zgrep '^#' $outfile; zgrep -v '^#' $outfile | sort -k1,1n -k2n ) | gzip -c >$outfile.1 && mv $outfile.1 $outfile`;
    } else {
        `(grep '^#' $outfile; grep -v '^#' $outfile | sort -k1,1n -k2n) >$outfile.1 && mv $outfile.1 $outfile`;
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
    
    foreach my $gt (@$gts) {
    
        # unphased or phased.
        $allele->{$ref} += 2 if ($gt eq '0/0' || $gt eq '0|0');
        $allele->{$alt} += 2 if ($gt eq '1/1' || $gt eq '1|1');
        if ($gt eq '0/1' || $gt eq '0|1' || $gt eq '1|0') {
            $allele->{$ref}++; 
            $allele->{$alt}++;
            $allele->{'het'}++;
        }
    }
    my $n    = ($allele->{$ref} + $allele->{$alt}) / 2;
    my $p    = $allele->{$ref} / (2 * $n);
    my $Ho   = $allele->{'het'} / $n;
    my $He   = 2 * $p * (1 - $p);
    my $Fis  = $He > 0 ? (1 - $Ho/$He) : 'NAN';
    my $flag = $p < 0.5 ? $ref : $alt;
    my $maf  = $p < 0.5 ? $p : (1 - $p);
    
    return($Ho, $He, $Fis, $maf, $flag);
    
}

sub filter {

    my $vcf = shift;
    foreach my $vcfline (@$vcf) {
    #################################################################
    
    #VCF format:
    #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Indiv_genotypes ...
    #Genotypes:GT:PL:DP:SP:GQ ...
    
    #Filter order: biallelic -> Q -> depth -> GQ -> coverage -> local Ho or Fis or MAF
    #-> missing rate -> global Ho or Fis or MAF.
    
    #################################################################
    next if $vcfline =~ /^#/; 
    chomp $vcfline;
    
    if (1) {
        lock($snp_total);
        $snp_total++; # Number of total SNPs.
    }
    
    my $tot_miss = 0; # Missing rate for a site.
    my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @genos) = split(/\t/, $vcfline);
    my ($fmt, @filtered, @gts_total);
    undef @filtered;
    undef @gts_total;
    my $cnt_Fis          = 0;
    my $cnt_Ho           = 0;
    my $cnt_gmaf         = 0;
    my $cnt_lmaf->{$ref} = 0;
       $cnt_lmaf->{$alt} = 0;
    my $mono             = 0;
    next if $alt =~ /,|\./;                         # Skip non-biallelic loci.
    next if (defined $minQ && $qual < $minQ);       # Minmum quality.
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
            if (not defined $geno[$fmt->{'DP'}]) { die "DP field is requred!\n";}
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
            
            my $DP   = $geno[$fmt->{'DP'}];
            my $GQ   = $geno[$fmt->{'GQ'}] if $fmt->{'GQ'};
            
            #Other filter...
            ###### Depth ######
            if ((defined $minDP && $DP < $minDP) or (defined $maxDP && $DP > $maxDP)) {
                
                $geno[$fmt->{'GT'}] = './.';
                $miss++;
                $tot_miss++;
                push @filtered, join(":", @geno);
                next;
            }
            
            ###### GQ ######
            if (defined $minGQ && $GQ < $minGQ) {
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
            last if $cov_ratio < $cov;
        }
        
        #next if $miss == $total;
        
        ###### Local ######
        if ((!$global && (defined $het or defined $fis)) or $l_maf) {
        
            my ($Ho, $He, $Fis, $maf, $flag) = calc_stat(\@gts, $ref, $alt); #Each population.
            $maf = ($l_maf && $l_maf > 1) ? 2*$l_N*$maf : $maf; # allele count or freq.
            ###### Het ######
            if (defined $het && $Ho > $het) {
                $cnt_Ho++;
                last;
            }
        
            ###### Fis ######
            if (defined $fis && abs($Fis) > $fis) {
                $cnt_Fis++;
                last;
            }
        
            ###### MAF ######
            $cnt_lmaf->{$flag}++ if $l_maf && $maf < $l_maf;
        }
        
    }
    
    next if ($recode  > 0 || $cnt_Ho  > 0 || $cnt_Fis > 0);
    
    ###### Global missing rate #######
    my $g_N = $tot_indiv - $tot_miss; # global num of individuals.
    if (defined $tot_cov) {
        my $cov_ratio = $tot_cov > 1 ? $g_N : $g_N/$tot_indiv;
        next if $cov_ratio < $tot_cov;
    }
    
    ###### Global ######
    if ($global or $g_maf) {
        
        my ($Ho, $He, $Fis, $maf, $flag) = calc_stat(\@gts_total, $ref, $alt);
        $maf = ($g_maf && $g_maf > 1) ? 2*$g_N*$maf : $maf;
        ###### Het ######
        if ($global && defined $het && $Ho > $het) {
            $cnt_Ho++;
            next;
        }
        
        ###### Fis ######
        if ($global && defined $fis && abs($Fis) > $fis) {
            $cnt_Fis++;
            next;
        }
        
        ###### MAF ######
        $cnt_gmaf++ if defined $g_maf && $maf < $g_maf;
        if (not defined $l_maf) {
            $cnt_lmaf->{$ref} = scalar(@order);
            $cnt_lmaf->{$alt} = scalar(@order);
        }
        
    }
       
    next if ($cnt_Ho   > 0 || $cnt_Fis  > 0);
    next if $cnt_gmaf > 0 && ($cnt_lmaf->{$ref} == scalar(@order) or $cnt_lmaf->{$alt} == scalar(@order));
    
    ## check if is snp.
    next if scalar(uniq @gts_total) == 1;
    
    ###### Print each snp sites ######
    lock(@out);
    lock($cnt);
    push @out, join("\t", $chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @filtered) . "\n";
    $cnt++; # SNP number.
    #printf STDERR "Filtering: %-1s [Retained loci: %d]\r", $x, $cnt;
    #printf STDERR "Filtering: retained %d of total %d loci, %d individuals of %d populations.\r", $cnt, $snp_total, $tot_indiv, $num_pops;
    }
}