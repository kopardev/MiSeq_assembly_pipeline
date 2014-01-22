#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Config::General;
use Data::Dumper;
use lib qw(/usr/global/blp/perllib);
use Vutil;

sub usage;
sub getReadStats;
sub getperc;


if (@ARGV==0) {
    usage();
}


my ($read1,$read2,$sampleName,$configFileName,$ncpus,$bt2IndexListFile,@bt2IndexList,$unaligned,$minReadLength,$nodust,$keepFiles,$help);
my ($windowsize,$window_quality_threshold,$avg_read_quality_threshold);


my $result = GetOptions ( "c=s" => \$configFileName,
			  "1=s" => \$read1,
			  "2=s" => \$read2,
			  "n=i" => \$ncpus,
                          "s=s" => \$sampleName,
			  "m=i" => \$minReadLength,
			  "nodust" => \$nodust,
		          "k" => \$keepFiles,
			  "w=i" => \$windowsize,
			  "wqt=i" => \$window_quality_threshold,
			  "arqt=i" => \$avg_read_quality_threshold,
			  "h|help" => \$help);


$configFileName="/home/vnkoparde/pipelines_dev/MiSeq_assembly_pipeline/config.txt" unless (defined $configFileName);
usage() if ($help);
Vutil::fileCheck($configFileName,"Check the config file!");
if (not defined $read1 or not defined $read2 or not defined $sampleName) {
    print "Read1 or Read2 OR sampleName not defined\n";
    exit;
}
Vutil::fileCheck($read1,"Check read1 file!");
Vutil::fileCheck($read2,"Check read2 file!");
$minReadLength=50 unless (defined $minReadLength);
$ncpus=1 unless (defined $ncpus);
$windowsize=9 unless (defined $windowsize);
$window_quality_threshold=25 unless (defined $window_quality_threshold);
$avg_read_quality_threshold=30 unless (defined $avg_read_quality_threshold);

my $cfg=new Config::General($configFileName);
my %cfgHash=$cfg->getall;

open(R,">${sampleName}.report.txt");
print R "#MiSeq Assembly Pipeline Report\n\n";
print R "Time=".localtime()."\n";
print R "Read1=$read1\n";
print R "Read2=$read1\n";

my %initialStatsR1=getReadStats($read1);
my %initialStatsR2=getReadStats($read2);

print R "\n#Initial Stats\tRead1\tRead2\n";
print R "Number_of_reads\t$initialStatsR1{Number_of_reads}\t$initialStatsR2{Number_of_reads}\n";
print R "Number_of_bases\t$initialStatsR1{Number_of_bases}\t$initialStatsR2{Number_of_bases}\n";
print R "Minimum_read_length\t$initialStatsR1{Minimum_read_length}\t$initialStatsR2{Minimum_read_length}\n";
print R "Maximum_read_length\t$initialStatsR1{Maximum_read_length}\t$initialStatsR2{Maximum_read_length}\n";
print R "Average_read_length\t$initialStatsR1{Average_read_length}\t$initialStatsR2{Average_read_length}\n";
print R "Average_read_quality\t$initialStatsR1{Average_read_quality}\t$initialStatsR2{Average_read_quality}\n";

#print Data::Dumper->Dumper(\%initialStatsR1),"\n";
#print Data::Dumper->Dumper(\%initialStatsR2),"\n";
my $qtread1=$sampleName."_qt_R1.fastq.gz";
my $qtread2=$sampleName."_qt_R2.fastq.gz";
my $cmd="$cfgHash{fastq_trim_by_qual_pe} $read1 $read2 $qtread1 $qtread2 $windowsize $window_quality_threshold $minReadLength $avg_read_quality_threshold";
print $cmd."\n";
system($cmd);
my %step1StatsR1=getReadStats($qtread1);
my %step1StatsR2=getReadStats($qtread2);

print R "\n#After Quality trimming\tRead1\tRead2\n";
print R "Number_of_reads\t$step1StatsR1{Number_of_reads} (".getperc($step1StatsR1{Number_of_reads},$initialStatsR1{Number_of_reads}).")\t$step1StatsR2{Number_of_reads} (".getperc($step1StatsR2{Number_of_reads},$initialStatsR2{Number_of_reads}).")\n";
print R "Number_of_bases\t$step1StatsR1{Number_of_bases} (".getperc($step1StatsR1{Number_of_bases},$initialStatsR1{Number_of_bases}).")\t$step1StatsR2{Number_of_bases} (".getperc($step1StatsR2{Number_of_bases},$initialStatsR2{Number_of_bases}).")\n";
print R "Minimum_read_length\t$step1StatsR1{Minimum_read_length} (".getperc($step1StatsR1{Minimum_read_length},$initialStatsR1{Minimum_read_length}).")\t$step1StatsR2{Minimum_read_length} (".getperc($step1StatsR2{Minimum_read_length},$initialStatsR2{Minimum_read_length}).")\n";
print R "Maximum_read_length\t$step1StatsR1{Maximum_read_length} (".getperc($step1StatsR1{Maximum_read_length},$initialStatsR1{Maximum_read_length}).")\t$step1StatsR2{Maximum_read_length} (".getperc($step1StatsR2{Maximum_read_length},$initialStatsR2{Maximum_read_length}).")\n";
print R "Average_read_length\t$step1StatsR1{Average_read_length} (".getperc($step1StatsR1{Average_read_length},$initialStatsR1{Average_read_length}).")\t$step1StatsR2{Average_read_length} (".getperc($step1StatsR2{Average_read_length},$initialStatsR2{Average_read_length}).")\n";
print R "Average_read_quality\t$step1StatsR1{Average_read_quality} (".getperc($step1StatsR1{Average_read_quality},$initialStatsR1{Average_read_quality}).")\t$step1StatsR2{Average_read_quality} (".getperc($step1StatsR2{Average_read_quality},$initialStatsR2{Average_read_quality}).")\n";

my $lcread1=$sampleName."_qtlc_R1.fastq.gz";
my $lcread2=$sampleName."_qtlc_R2.fastq.gz";
my $lcout=$sampleName.".lc.fastq";
$cmd="$cfgHash{sga} preprocess -m $minReadLength --dust -p 1 $qtread1 $qtread2 > $lcout";
print $cmd."\n";
system($cmd);
$cmd="$cfgHash{fastq_deinterleave} $lcout ${sampleName}_qtlc";
print $cmd."\n";
system($cmd);

my %step2StatsR1=getReadStats($lcread1);
my %step2StatsR2=getReadStats($lcread2);

print R "\n#After Quality trimming\tRead1\tRead2\n";
print R "Number_of_reads\t$step2StatsR1{Number_of_reads} (".getperc($step2StatsR1{Number_of_reads},$initialStatsR1{Number_of_reads}).")\t$step2StatsR2{Number_of_reads} (".getperc($step2StatsR2{Number_of_reads},$initialStatsR2{Number_of_reads}).")\n";
print R "Number_of_bases\t$step2StatsR1{Number_of_bases} (".getperc($step2StatsR1{Number_of_bases},$initialStatsR1{Number_of_bases}).")\t$step2StatsR2{Number_of_bases} (".getperc($step2StatsR2{Number_of_bases},$initialStatsR2{Number_of_bases}).")\n";
print R "Minimum_read_length\t$step2StatsR1{Minimum_read_length} (".getperc($step2StatsR1{Minimum_read_length},$initialStatsR1{Minimum_read_length}).")\t$step2StatsR2{Minimum_read_length} (".getperc($step2StatsR2{Minimum_read_length},$initialStatsR2{Minimum_read_length}).")\n";
print R "Maximum_read_length\t$step2StatsR1{Maximum_read_length} (".getperc($step2StatsR1{Maximum_read_length},$initialStatsR1{Maximum_read_length}).")\t$step2StatsR2{Maximum_read_length} (".getperc($step2StatsR2{Maximum_read_length},$initialStatsR2{Maximum_read_length}).")\n";
print R "Average_read_length\t$step2StatsR1{Average_read_length} (".getperc($step2StatsR1{Average_read_length},$initialStatsR1{Average_read_length}).")\t$step2StatsR2{Average_read_length} (".getperc($step2StatsR2{Average_read_length},$initialStatsR2{Average_read_length}).")\n";
print R "Average_read_quality\t$step2StatsR1{Average_read_quality} (".getperc($step2StatsR1{Average_read_quality},$initialStatsR1{Average_read_quality}).")\t$step2StatsR2{Average_read_quality} (".getperc($step2StatsR2{Average_read_quality},$initialStatsR2{Average_read_quality}).")\n";

close R;
exit;


print $qtread1," ",$qtread2,"\n";exit;
#my $nreads_lq=0;
#my $nreads_lc=0;
#
#if ( $dust==1 ) {
#    my $HQread1=getHQname($read1);
#    my $HQread2=getHQname($read2);
#    my $cmd="export QC_PRINTHQ=1 && /home/vnkoparde/bin/fastq_getQCStatsPE $read1 $read2";
#    system($cmd) unless ( -e $HQread1 and -e $HQread2 );
#    $nreads_lq=$nreads-getnreads($HQread1);
#    print "/home/vnkoparde/opt/sga/bin/sga -m $minReadLength --dust -p 1 $HQread1 $HQread2\n";
#    my $cmd1="/home/vnkoparde/opt/sga/bin/sga dust -m $minReadLength --dust -p 1 $HQread1 $HQread2 |gzip -c - >  ${sampleName}_lqlc_filtered.interleaved.fastq.gz";
#    my $cmd2="/home/vnkoparde/bin/fastq_deinterleave ${sampleName}_lqlc_filtered.interleaved.fastq.gz ${sampleName}_lqlc_filtered";
#    system($cmd1) unless ( -e "${sampleName}_lqlc_filtered.interleaved.fastq.gz");
#    system($cmd2) unless ( -e "${sampleName}_lqlc_filtered_R1.fastq.gz" and -e "${sampleName}_lqlc_filtered_R2.fastq.gz");
#    $read1="${sampleName}_lqlc_filtered_R1.fastq.gz";
#    $read2="${sampleName}_lqlc_filtered_R2.fastq.gz";
#    $nreads_lc=$nreads-$nreads_lq-getnreads($read1);
#}
#
#open BT2ILF, "<$bt2IndexListFile";
#while (my $line=<BT2ILF>) {
#    chomp $line;
#    push @bt2IndexList,$line;
#}
#close BT2ILF;
#
#my $cmd;
#for my $bt2Index (@bt2IndexList) {
#    my @tmp=split/\//,$bt2Index;
#    my $indexName=$tmp[-1];
#    my $gi2tax=$bt2Index.".gi2tax.gz";
#    $cmd="$bowtie2Bin $bowtie2parameters -p $ ncpus -x $bt2Index -1 $read1 -2 $read2 | $samtoolsBin view -f2 -bS - > ${sampleName}.${indexName}.bam";
#    system($cmd) unless ( -e "${sampleName}.${indexName}.bam" );
#    $cmd="perl /home/vnkoparde/scripts/bam2percaln.v2.pl ${sampleName}.${indexName}.bam | perl /home/vnkoparde/scripts/percalnAddTaxa.pl - $gi2tax | /bin/gzip -c - > ${sampleName}.${indexName}.percaln.tax.gz";
#    system($cmd) unless ( -e "${sampleName}.${indexName}.percaln.tax.gz" );
#}
#$cmd="zcat ${sampleName}.*.percaln.tax.gz | sort -k1,1 -k6,6gr | perl /home/vnkoparde/scripts/percaln_remdup.pl - | perl /home/vnkoparde/scripts/percalntax2profile.pl - $nreads $nreads_lq $nreads_lc $sampleName";
#system($cmd);

sub usage {
print <<EOF;
MiSeq Assembly Pipeline

Note:
1. Only works with paired end Illumina data.

List of output files:

Author: Vishal N. Koparde, Ph. D.
Created: 140122
Modified: 140122

options:
--1 read1.fastq or read1.fastq.gz
--2 read2.fastq or read2.fastq.gz
--s name of the sample
--n number of cpus
--nodust do not perform dusting (low complexity filtering)
--k keep intermediate files
--m minimum read length to consider (default=50)
--w sliding window width for quality trimming (default=9)
--wqt sliding window quality threshold (default=25)
--arqt average read quality threshold (default=30)

EOF
exit 1;
}

sub getReadStats {
    my ($fq)=@_;
    my $fq_stats="${fq}.numreads+";
    my $nreads=`$cfgHash{"fastq_num_reads+"} $fq > $fq_stats`;
    my $tmp=Config::General->new($fq_stats);
    return $tmp->getall;    
}

sub getperc {
    my ($num,$den)=@_;
    my $p=100.0*$num/$den;
    my $q=sprintf "%.2f",$p;
    return $q;
}
