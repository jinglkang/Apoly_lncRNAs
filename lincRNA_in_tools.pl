#!/usr/bin/perl -w
use strict;
use warnings;

# CPC_lncRNA_gene.txt; CPAT_lncRNA_gene.txt; FEELnc_lncRNA_gene.txt
# final lincRNAs file: Candidate_lncRNA_genes.classifier.final.filter_reads.txt
my %CPC   = &build_hash("CPC_lncRNA_gene.txt");
my %CPAT  = &build_hash("CPAT_lncRNA_gene.txt");
my %FEELnc= &build_hash("FEELnc_lncRNA_gene.txt");

my $lincRNA="classifier/Candidate_lncRNA_genes.classifier.final.filter_reads.txt";
open LINC, $lincRNA or die "can not open $lincRNA\n";
while (<LINC>) {
        chomp;
        my @a=split /\t/;
        next if /^LncRNA_gene/;
        my ($a1, $a2, $a3);
        if ($a[4] eq "intergenic") {
                $CPC{$a[0]}?($a1="CPC"):($a1="-");
                $CPAT{$a[0]}?($a2="CPAT"):($a2="-");
                $FEELnc{$a[0]}?($a3="FEELnc"):($a3="-");
                my $value=$a1.";".$a2.";".$a3;
                print "$a[0]\tintergenic\t$value\n";
        }
}

sub build_hash {
        my ($file)=@_;
        my %hash;
        open FIL, $file or die "can not open $file\n";
        while (<FIL>) {
                chomp;
                $hash{$_}++;
        }
        return %hash;
}
