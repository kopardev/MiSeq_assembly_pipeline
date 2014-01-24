
#! /usr/bin/perl -W

use strict;
use Bio::SeqIO;
use lib qw(/usr/global/blp/perllib);
use Vutil;

my $fastaFile=shift;

if(!defined($fastaFile)) { 
    die "USAGE: perl $0 fastafile\n";
}

Vutil::fileCheck($fastaFile,"Cannot calculate genome statistics!");

my $in=Bio::SeqIO->new(-file=>$fastaFile,-format=>'fasta');
my %seqGC;
my %seqN;
my %seqLength;
my $genomeSize;
my $largestContigSize;
my $smallestContigSize;
my $ncontigs;
my $nContigs10kPlus;
my $sizeContigs10kPlus;
my $sumGC;
my $sumN;

$ncontigs=0;
$smallestContigSize=10000000000000;
$largestContigSize=0;
$genomeSize=0;
$nContigs10kPlus=0;
$sizeContigs10kPlus=0;
$sumGC=0;
$sumN=0;
while (my $s=$in->next_seq) {
    $ncontigs+=1;
    my @tmp=($s->seq)=~/[GC]/ig;
    my $gc=scalar @tmp;
    @tmp=($s->seq)=~/[N]/ig;
    my $n=scalar @tmp;
    $seqGC{$s->id}=$gc;
    $seqN{$s->id}=$n;
    $sumGC+=$gc;
    $sumN+=$n;
    $seqLength{$s->id}=$s->length;
    $genomeSize+=$s->length;
    if ($s->length>10000) {
	$nContigs10kPlus+=1;
	$sizeContigs10kPlus+=$s->length;
    }
    $smallestContigSize=($smallestContigSize>$s->length)?$s->length:$smallestContigSize;
    $largestContigSize=($largestContigSize<$s->length)?$s->length:$largestContigSize;
}
exit if $genomeSize==0;
my $avgContigSize=int($genomeSize/$ncontigs);
my $percentGC=sprintf "%.2f",$sumGC*100.0/$genomeSize;
my $percentN=sprintf "%.2f",$sumN*100.0/$genomeSize;
my $percentnContigs10kPlus=sprintf "%.2f",$nContigs10kPlus*100.0/$ncontigs;
my $percentsizeContigs10kPlus=sprintf "%.2f",$sizeContigs10kPlus*100.0/$genomeSize;

my ($n50Size,$n50Number,$n50AvgContigSize,$n50Found);
my ($n90Size,$n90Number,$n90AvgContigSize,$n90Found);

$n50Size=0;
$n50Number=0;
$n50AvgContigSize=0;
$n90Size=0;
$n90Number=0;
$n90AvgContigSize=0;

my $n50Target;
my $n90Target;
$n50Target=0.5*$genomeSize;
$n90Target=0.9*$genomeSize;

$n50Found=0;
$n90Found=0;
my $dummySize;
my $dummynContigs;
$dummySize=0;
$dummynContigs=0;

foreach my $sid (sort { $seqLength{$b} <=> $seqLength{$a} } keys %seqLength ) {
    last if ($n90Found!=0 and $n50Found != 0);
    $dummySize+=$seqLength{$sid};
    $dummynContigs+=1;
    #print "$dummynContigs\t$dummySize\t$n50Target\t$n90Target\t$n50Found\t$n90Found\n";
    if ($n50Found==0) {
	if ($dummySize>$n50Target) {
	    $n50Found=1;
	    $n50Size=$seqLength{$sid};
	    $n50Number=$dummynContigs;
	    $n50AvgContigSize=int($dummySize/$dummynContigs);
	}
    }
    if ($n90Found==0) {
	if ($dummySize>$n90Target) {
	    $n90Found=1;
	    $n90Size=$seqLength{$sid};
	    $n90Number=$dummynContigs;
	    $n90AvgContigSize=int($dummySize/$dummynContigs);
	}
    }
}

#print "$ncontigs\t$largestContigSize\t$smallestContigSize\t$avgContigSize\t$n50Size\t$n50Number\t$n50AvgContigSize\t$n90Size\t$n90Number\t$n90AvgContigSize\t$genomeSize\t$percentGC\t$sumN\t$percentN\t$nContigs10kPlus\t$percentnContigs10kPlus\t$sizeContigs10kPlus\t$percentsizeContigs10kPlus\n";
print "Number_of_contigs\t$ncontigs\n";
print "Largest_contig_size\t$largestContigSize\n";
print "Smallest_contig_size\t$smallestContigSize\n";
print "Average_contig_size\t$avgContigSize\n";
print "N50Size\t$n50Size\n";
print "N50Number\t$n50Number\n";
print "N50Avg_contig_size\t$n50AvgContigSize\n";
print "N90Size\t$n90Size\n";
print "N90Number\t$n90Number\n";
print "N90Avg_contig_size\t$n90AvgContigSize\n";
print "GenomeSize\t$genomeSize\n";
print "GC\t$percentGC\n";
print "Number_of_Ns\t$sumN\n";
print "Percent_Ns\t$percentN\n";
print "Number_of_contigs>10k\t$nContigs10kPlus\n";
print "Percent_of_contigs>10k\t$percentnContigs10kPlus\n";
print "Total_size_of_contigs>10k\t$sizeContigs10kPlus\n";
print "Percent_of_genome_in_contigs>10k\t$percentsizeContigs10kPlus\n";
