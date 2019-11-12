#!/usr/bin/perl -w

# A perl script to build config file for SRS Bam Merge
# Contact: Shihcheng.Guo@Gmail.com
# Version 1.3
# Go to http://sra.dnanexus.com/studies/SRP028600/samples
# 11/12/2019
# Select SRS and Click Related RUNS then get the Table as the input

use strict;
use warnings;
use Cwd;

my $file=shift @ARGV;
my %SRR;
my %SRX;
open F,$file;
while(<F>){
chomp;
my @line=split/\s+/;
push @{$SRX{$line[2]}}, $line[0];
print "$line[0]\t$line[2]\n";
}

foreach my $srx(sort keys %SRX){
my @SRR=@{$SRX{$srx}};
my %beta;
foreach my $srr(@SRR){
print "$srr\n";
open F,"$srr\_1_val_1_bismark_bt2_pe.nonCG_filtered.bismark.cov" || die "cannot find $srr";
while(<F>){
my @line=split/\s+/;
my $pos="$line[0]:$line[1]-$line[2]";
$beta{$pos}{"M"}+=$line[4];
$beta{$pos}{"U"}+=$line[5];
}
}
open OUT,">$srx.bedgraph";
foreach my $pos(sort keys %beta){
my $beta=$beta{$pos}{"M"}/($beta{$pos}{"M"}+$beta{$pos}{"U"});
my ($CHR,$START,$END)=split/:|-/,$pos;
print OUT "$CHR\t$START\t$END\t$beta\t$beta{$pos}{'M'}\t$beta{$pos}{'U'}\n";
}
close OUT;
}
