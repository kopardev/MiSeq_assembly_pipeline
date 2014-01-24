#!/usr/bin/perl

#Vishal Koparde
#v1

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
my ($windowsize,$window_quality_threshold,$avg_read_quality_threshold,$noqt,$assemblyFile,$nonodup,$nopat,$nocov);


my $result = GetOptions ( "c=s" => \$configFileName,
			  "1=s" => \$read1,
			  "2=s" => \$read2,
			  "n=i" => \$ncpus,
                          "s=s" => \$sampleName,
			  "m=i" => \$minReadLength,
			  "nonodup" => \$nonodup,
			  "noqt" => \$noqt,
			  "nodust" => \$nodust,
			  "nopat" => \$nopat,
			  "nocov" => \$nocov,
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
$assemblyFile=$sampleName."_contigs.fasta";

my $cfg=new Config::General($configFileName);
my %cfgHash=$cfg->getall;
my @todelete;
my $cmd;

open(R,">${sampleName}.report.txt");
print R "#MiSeq Assembly Pipeline Report\n\n";
print R "Time=".localtime()."\n";
print R "Read1=$read1\n";
print R "Read2=$read1\n";
my $dummy;
$dummy=(defined $nonodup)?"FALSE":"TRUE";
print R "RemoveDuplicates=$dummy\n";
$dummy=(defined $noqt)?"FALSE":"TRUE";
print R "QualityTrimming=$dummy\n";
$dummy=(defined $nodust)?"FALSE":"TRUE";
print R "LowComplexityFiltering=$dummy\n";
$dummy=(defined $nopat)?"FALSE":"TRUE";
print R "PolyA/T_and_NFiltering=$dummy\n";

my %initialStatsR1=getReadStats($read1);
my %initialStatsR2=getReadStats($read2);

#print R "\n#Initial Stats\tRead1\tRead2\n";
#print R "Number_of_reads\t$initialStatsR1{Number_of_reads}\t$initialStatsR2{Number_of_reads}\n";
#print R "Number_of_bases\t$initialStatsR1{Number_of_bases}\t$initialStatsR2{Number_of_bases}\n";
#print R "Minimum_read_length\t$initialStatsR1{Minimum_read_length}\t$initialStatsR2{Minimum_read_length}\n";
#print R "Maximum_read_length\t$initialStatsR1{Maximum_read_length}\t$initialStatsR2{Maximum_read_length}\n";
#print R "Average_read_length\t$initialStatsR1{Average_read_length}\t$initialStatsR2{Average_read_length}\n";
#print R "Average_read_quality\t$initialStatsR1{Average_read_quality}\t$initialStatsR2{Average_read_quality}\n";

#print Data::Dumper->Dumper(\%initialStatsR1),"\n";
#print Data::Dumper->Dumper(\%initialStatsR2),"\n";

# remove possible PCR duplicates

my $nodupRead1;
my $nodupRead2;

if (defined $nonodup) {
    $nodupRead1=$read1;
    $nodupRead2=$read2;
} else {
    my @tmp1=split/\./,$read1;
    my @tmp2=split/\./,$read2;
    pop @tmp1 if $tmp1[-1] eq "gz";
    pop @tmp2 if $tmp2[-1] eq "gz";
    pop @tmp1 if $tmp1[-1] eq "fastq";
    pop @tmp2 if $tmp2[-1] eq "fastq";
    push @tmp1,("nodup","fastq","gz");
    push @tmp2,("nodup","fastq","gz");
    $nodupRead1=join(".",@tmp1);
    $nodupRead2=join(".",@tmp2);
    $ENV{QC_PRINTNODUP}=1;
    $cmd="$cfgHash{fastq_getQCStatsPE} $read1 $read2";
    print $cmd."\n";
    system($cmd);
    push @todelete,$nodupRead1;
    push @todelete,$nodupRead2;
}

my %step1StatsR1=getReadStats($nodupRead1);
my %step1StatsR2=getReadStats($nodupRead2);

#print R "\n#After Duplicate removal\tRead1\tRead2\n";
#print R "Number_of_reads\t$step1StatsR1{Number_of_reads} (".getperc($step1StatsR1{Number_of_reads},$initialStatsR1{Number_of_reads})."\%)\t$step1StatsR2{Number_of_reads} (".getperc($step1StatsR2{Number_of_reads},$initialStatsR2{Number_of_reads})."\%)\n";
#print R "Number_of_bases\t$step1StatsR1{Number_of_bases} (".getperc($step1StatsR1{Number_of_bases},$initialStatsR1{Number_of_bases})."\%)\t$step1StatsR2{Number_of_bases} (".getperc($step1StatsR2{Number_of_bases},$initialStatsR2{Number_of_bases})."\%)\n";
#print R "Minimum_read_length\t$step1StatsR1{Minimum_read_length}\t$step1StatsR2{Minimum_read_length}\n";
#print R "Maximum_read_length\t$step1StatsR1{Maximum_read_length}\t$step1StatsR2{Maximum_read_length}\n";
#print R "Average_read_length\t$step1StatsR1{Average_read_length}\t$step1StatsR2{Average_read_length}\n";
#print R "Average_read_quality\t$step1StatsR1{Average_read_quality}\t$step1StatsR2{Average_read_quality}\n";

# do quality trimming from 3'

my $qtread1;
my $qtread2;

if (defined $noqt) {
    $qtread1=$nodupRead1;
    $qtread2=$nodupRead2;
} else {
    $qtread1=$sampleName."_qt_R1.fastq.gz";
    $qtread2=$sampleName."_qt_R2.fastq.gz";
    $cmd="$cfgHash{fastq_trim_by_qual_pe} $nodupRead1 $nodupRead2 $qtread1 $qtread2 $windowsize $window_quality_threshold $minReadLength $avg_read_quality_threshold";
    print $cmd."\n";
    system($cmd);
    push @todelete,$qtread1;
    push @todelete,$qtread2;
}

my %step2StatsR1=getReadStats($qtread1);
my %step2StatsR2=getReadStats($qtread2);

#print R "\n#After Quality trimming\tRead1\tRead2\n";
#print R "Number_of_reads\t$step2StatsR1{Number_of_reads} (".getperc($step2StatsR1{Number_of_reads},$initialStatsR1{Number_of_reads})."\%)\t$step2StatsR2{Number_of_reads} (".getperc($step2StatsR2{Number_of_reads},$initialStatsR2{Number_of_reads})."\%)\n";
#print R "Number_of_bases\t$step2StatsR1{Number_of_bases} (".getperc($step2StatsR1{Number_of_bases},$initialStatsR1{Number_of_bases})."\%)\t$step2StatsR2{Number_of_bases} (".getperc($step2StatsR2{Number_of_bases},$initialStatsR2{Number_of_bases})."\%)\n";
#print R "Minimum_read_length\t$step2StatsR1{Minimum_read_length}\t$step2StatsR2{Minimum_read_length}\n";
#print R "Maximum_read_length\t$step2StatsR1{Maximum_read_length}\t$step2StatsR2{Maximum_read_length}\n";
#print R "Average_read_length\t$step2StatsR1{Average_read_length}\t$step2StatsR2{Average_read_length}\n";
#print R "Average_read_quality\t$step2StatsR1{Average_read_quality}\t$step2StatsR2{Average_read_quality}\n";

# remove low complexity using sga

my $lcread1;
my $lcread2;
if (defined $nodust) {
    $lcread1=$qtread1;
    $lcread2=$qtread2;
} else {
    $lcread1=$sampleName."_lc_R1.fastq.gz";
    $lcread2=$sampleName."_lc_R2.fastq.gz";
    my $lcout=$sampleName.".lc.fastq";
    $cmd="$cfgHash{sga} preprocess -m $minReadLength --dust -p 1 $qtread1 $qtread2 > $lcout";
    print $cmd."\n";
    system($cmd);
    $cmd="$cfgHash{fastq_deinterleave} $lcout ${sampleName}_lc";
    print $cmd."\n";
    system($cmd);
    push @todelete,$lcread1;
    push @todelete,$lcread2;
    push @todelete,$lcout;
}

my %step3StatsR1=getReadStats($lcread1);
my %step3StatsR2=getReadStats($lcread2);

#print R "\n#After Low Complexity filtering\tRead1\tRead2\n";
#print R "Number_of_reads\t$step3StatsR1{Number_of_reads} (".getperc($step3StatsR1{Number_of_reads},$initialStatsR1{Number_of_reads})."\%)\t$step3StatsR2{Number_of_reads} (".getperc($step3StatsR2{Number_of_reads},$initialStatsR2{Number_of_reads})."\%)\n";
#print R "Number_of_bases\t$step3StatsR1{Number_of_bases} (".getperc($step3StatsR1{Number_of_bases},$initialStatsR1{Number_of_bases})."\%)\t$step3StatsR2{Number_of_bases} (".getperc($step3StatsR2{Number_of_bases},$initialStatsR2{Number_of_bases})."\%)\n";
#print R "Minimum_read_length\t$step3StatsR1{Minimum_read_length}\t$step3StatsR2{Minimum_read_length}\n";
#print R "Maximum_read_length\t$step3StatsR1{Maximum_read_length}\t$step3StatsR2{Maximum_read_length}\n";
#print R "Average_read_length\t$step3StatsR1{Average_read_length}\t$step3StatsR2{Average_read_length}\n";
#print R "Average_read_quality\t$step3StatsR1{Average_read_quality}\t$step3StatsR2{Average_read_quality}\n";


# remove polyA/T from 5' and 3' ends and reads with too many Ns using fqtrim

my $patread1;
my $patread2;
if (defined $nopat) {
    $patread1=$lcread1;
    $patread2=$lcread2;
} else {
    my @tmp1=split/\./,$lcread1;
    my @tmp2=split/\./,$lcread2;
    pop @tmp1 if $tmp1[-1] eq "gz";
    pop @tmp2 if $tmp2[-1] eq "gz";
    pop @tmp1 if $tmp1[-1] eq "fastq";
    pop @tmp2 if $tmp2[-1] eq "fastq";
    push @tmp1,("pat","fastq","gz");
    push @tmp2,("pat","fastq","gz");
    $patread1=join(".",@tmp1);
    $patread2=join(".",@tmp2);
    $cmd="$cfgHash{fqtrim} -p $ncpus -B -P33 -l $minReadLength -o pat.fastq.gz $lcread1,$lcread2";
    print $cmd."\n";
    system($cmd);
    push @todelete,$patread1;
    push @todelete,$patread2;
}

my %step4StatsR1=getReadStats($patread1);
my %step4StatsR2=getReadStats($patread2);

#print R "\n#After Poly A/T trimming and N-filtering\tRead1\tRead2\n";
#print R "Number_of_reads\t$step4StatsR1{Number_of_reads} (".getperc($step4StatsR1{Number_of_reads},$initialStatsR1{Number_of_reads})."\%)\t$step4StatsR2{Number_of_reads} (".getperc($step4StatsR2{Number_of_reads},$initialStatsR2{Number_of_reads})."\%)\n";
#print R "Number_of_bases\t$step4StatsR1{Number_of_bases} (".getperc($step4StatsR1{Number_of_bases},$initialStatsR1{Number_of_bases})."\%)\t$step4StatsR2{Number_of_bases} (".getperc($step4StatsR2{Number_of_bases},$initialStatsR2{Number_of_bases})."\%)\n";
#print R "Minimum_read_length\t$step4StatsR1{Minimum_read_length}\t$step4StatsR2{Minimum_read_length}\n";
#print R "Maximum_read_length\t$step4StatsR1{Maximum_read_length}\t$step4StatsR2{Maximum_read_length}\n";
#print R "Average_read_length\t$step4StatsR1{Average_read_length}\t$step4StatsR2{Average_read_length}\n";
#print R "Average_read_quality\t$step4StatsR1{Average_read_quality}\t$step4StatsR2{Average_read_quality}\n";


print R "\n#Initial Stats\t\t\t";
print R "After Duplicate removal\t\t\t\t";
print R "After Quality trimming\t\t\t\t";
print R "After Low Complexity filtering\t\t\t\t";
print R "After Poly A/T trimming and N-filtering\t\t\t\n";

print R "#\tRead1\tRead2";
print R "\tRead1\t\tRead2\t";
print R "\tRead1\t\tRead2\t";
print R "\tRead1\t\tRead2\t";
print R "\tRead1\t\tRead2\n";

print R "Number_of_reads\t$initialStatsR1{Number_of_reads}\t$initialStatsR2{Number_of_reads}\t";
print R "$step1StatsR1{Number_of_reads}\t".getperc($step1StatsR1{Number_of_reads},$initialStatsR1{Number_of_reads})."\%\t$step1StatsR2{Number_of_reads}\t".getperc($step1StatsR2{Number_of_reads},$initialStatsR2{Number_of_reads})."\%\t";
print R "$step2StatsR1{Number_of_reads}\t".getperc($step2StatsR1{Number_of_reads},$initialStatsR1{Number_of_reads})."\%\t$step2StatsR2{Number_of_reads}\t".getperc($step2StatsR2{Number_of_reads},$initialStatsR2{Number_of_reads})."\%\t";
print R "$step3StatsR1{Number_of_reads}\t".getperc($step3StatsR1{Number_of_reads},$initialStatsR1{Number_of_reads})."\%\t$step3StatsR2{Number_of_reads}\t".getperc($step3StatsR2{Number_of_reads},$initialStatsR2{Number_of_reads})."\%\t";
print R "$step4StatsR1{Number_of_reads}\t".getperc($step4StatsR1{Number_of_reads},$initialStatsR1{Number_of_reads})."\%\t$step4StatsR2{Number_of_reads}\t".getperc($step4StatsR2{Number_of_reads},$initialStatsR2{Number_of_reads})."\%\n";

print R "Number_of_bases\t$initialStatsR1{Number_of_bases}\t$initialStatsR2{Number_of_bases}\t";
print R "$step1StatsR1{Number_of_bases}\t".getperc($step1StatsR1{Number_of_bases},$initialStatsR1{Number_of_bases})."\%\t$step1StatsR2{Number_of_bases}\t".getperc($step1StatsR2{Number_of_bases},$initialStatsR2{Number_of_bases})."\%\t";
print R "$step2StatsR1{Number_of_bases}\t".getperc($step2StatsR1{Number_of_bases},$initialStatsR1{Number_of_bases})."\%\t$step2StatsR2{Number_of_bases}\t".getperc($step2StatsR2{Number_of_bases},$initialStatsR2{Number_of_bases})."\%\t";
print R "$step3StatsR1{Number_of_bases}\t".getperc($step3StatsR1{Number_of_bases},$initialStatsR1{Number_of_bases})."\%\t$step3StatsR2{Number_of_bases}\t".getperc($step3StatsR2{Number_of_bases},$initialStatsR2{Number_of_bases})."\%\t";
print R "$step4StatsR1{Number_of_bases}\t".getperc($step4StatsR1{Number_of_bases},$initialStatsR1{Number_of_bases})."\%\t$step4StatsR2{Number_of_bases}\t".getperc($step4StatsR2{Number_of_bases},$initialStatsR2{Number_of_bases})."\%\n";

print R "Minimum_read_length\t$initialStatsR1{Minimum_read_length}\t$initialStatsR2{Minimum_read_length}\t";
print R "$step1StatsR1{Minimum_read_length}\t\t$step1StatsR2{Minimum_read_length}\t\t";
print R "$step2StatsR1{Minimum_read_length}\t\t$step2StatsR2{Minimum_read_length}\t\t";
print R "$step3StatsR1{Minimum_read_length}\t\t$step3StatsR2{Minimum_read_length}\t\t";
print R "$step4StatsR1{Minimum_read_length}\t\t$step4StatsR2{Minimum_read_length}\t\n";

print R "Maximum_read_length\t$initialStatsR1{Maximum_read_length}\t$initialStatsR2{Maximum_read_length}\t";
print R "$step1StatsR1{Maximum_read_length}\t\t$step1StatsR2{Maximum_read_length}\t\t";
print R "$step2StatsR1{Maximum_read_length}\t\t$step2StatsR2{Maximum_read_length}\t\t";
print R "$step3StatsR1{Maximum_read_length}\t\t$step3StatsR2{Maximum_read_length}\t\t";
print R "$step4StatsR1{Maximum_read_length}\t\t$step4StatsR2{Maximum_read_length}\t\n";

print R "Average_read_length\t$initialStatsR1{Average_read_length}\t$initialStatsR2{Average_read_length}\t";
print R "$step1StatsR1{Average_read_length}\t\t$step1StatsR2{Average_read_length}\t\t";
print R "$step2StatsR1{Average_read_length}\t\t$step2StatsR2{Average_read_length}\t\t";
print R "$step3StatsR1{Average_read_length}\t\t$step3StatsR2{Average_read_length}\t\t";
print R "$step4StatsR1{Average_read_length}\t\t$step4StatsR2{Average_read_length}\t\n";

print R "Average_read_quality\t$initialStatsR1{Average_read_quality}\t$initialStatsR2{Average_read_quality}\t";
print R "$step1StatsR1{Average_read_quality}\t\t$step1StatsR2{Average_read_quality}\t\t";
print R "$step2StatsR1{Average_read_quality}\t\t$step2StatsR2{Average_read_quality}\t\t";
print R "$step3StatsR1{Average_read_quality}\t\t$step3StatsR2{Average_read_quality}\t\t";
print R "$step4StatsR1{Average_read_quality}\t\t$step4StatsR2{Average_read_quality}\t\n";

print R "AssemblyRead1=$patread1\n";
print R "AssemblyRead2=$patread2\n";

# do not delete files used for assembly
@todelete = grep { $_ ne $patread1 } @todelete;
@todelete = grep { $_ ne $patread2 } @todelete;

# do the assembly

$cmd="$cfgHash{clc_assembler} -o $assemblyFile -p fb ss 100 400 -q -i $patread1 $patread2 -b 100 --cpus $ncpus -v";
system($cmd);

my $path=Vutil::getFileAbsolutePath(__FILE__);
Vutil::fileCheck("${path}/calculateAssemblyStats.pl","This script is required to calculate assembly stats!");
$cmd="perl ${path}/calculateAssemblyStats.pl $assemblyFile > ${sampleName}.assemblystats";
push @todelete,"${sampleName}.assemblystats";
print "$cmd\n";
system($cmd);
Vutil::fileCheck("${sampleName}.assemblystats","No assembly stats calculated!");
my $assemblyStatsTmp=Config::General->new("${sampleName}.assemblystats");
my %assemblyStats=$assemblyStatsTmp->getall;

print R "\n#Assembly Statistics\n";
print R "Number_of_contigs\t".$assemblyStats{"Number_of_contigs"}."\n";
print R "Largest_contig_size\t".$assemblyStats{"Largest_contig_size"}."\n";
print R "Smallest_contig_size\t".$assemblyStats{"Smallest_contig_size"}."\n";
print R "Average_contig_size\t".$assemblyStats{"Average_contig_size"}."\n";
print R "N50Size\t".$assemblyStats{"N50Size"}."\n";
print R "N50Number\t".$assemblyStats{"N50Number"}."\n";
print R "N50Avg_contig_size\t".$assemblyStats{"N50Avg_contig_size"}."\n";
print R "N90Size\t".$assemblyStats{"N90Size"}."\n";
print R "N90Number\t".$assemblyStats{"N90Number"}."\n";
print R "N90Avg_contig_size\t".$assemblyStats{"N90Avg_contig_size"}."\n";
print R "GenomeSize\t".$assemblyStats{"GenomeSize"}."\n";
print R "GC\t".$assemblyStats{"GC"}."\n";
print R "Number_of_Ns\t".$assemblyStats{"Number_of_Ns"}."\n";
print R "Percent_Ns\t".$assemblyStats{"Percent_Ns"}."\n";
print R "Number_of_contigs>10k\t".$assemblyStats{"Number_of_contigs>10k"}."\n";
print R "Percent_of_contigs>10k\t".$assemblyStats{"Percent_of_contigs>10k"}."\n";
print R "Total_size_of_contigs>10k\t".$assemblyStats{"Total_size_of_contigs>10k"}."\n";
print R "Percent_of_genome_in_contigs>10k\t".$assemblyStats{"Percent_of_genome_in_contigs>10k"}."\n";

my ($coverage,$nreadsAssembled,$nbasesAssembled,$percentReadsAssembled,$percentBasesAssembled);
$coverage=-1;
$nreadsAssembled=-1;
$nbasesAssembled=-1;
$percentReadsAssembled=-1;
$percentBasesAssembled=-1;

unless (defined $nocov) {
    $cmd="$cfgHash{clc_mapper} -o ${sampleName}.cas -q -p fb ss 100 400 -i $patread1 $patread2 -d $assemblyFile --cpus $ncpus";
    push @todelete,"${sampleName}.cas";
    print "$cmd\n";
    system($cmd);
    $cmd="$cfgHash{cas2sortedbam2} -i ${sampleName}.cas -u";
    push @todelete,"${sampleName}.bam";
    print "$cmd\n";
    system($cmd);
    $cmd="$cfgHash{bamGetBasicStats} ${sampleName}.bam > ${sampleName}.bam.basicstats";
    push @todelete,"${sampleName}.bam.basicstats";
    print "$cmd\n";
    system($cmd);
    my $bamStatsTmp=Config::General->new("${sampleName}.bam.basicstats");
    my %bamStats=$bamStatsTmp->getall;
    $coverage=sprintf "%.2f",$bamStats{BasesInMappedReads}/$assemblyStats{"GenomeSize"};
    $nreadsAssembled=int($bamStats{MappedReads});
    $percentReadsAssembled=sprintf "%.2f",100*$bamStats{MappedReads}/($step4StatsR1{Number_of_reads}+$step4StatsR2{Number_of_reads});
    $nbasesAssembled=int($bamStats{BasesInMappedReads});
    $percentBasesAssembled=sprintf "%.2f",100*$bamStats{BasesInMappedReads}/($step4StatsR1{Number_of_bases}+$step4StatsR2{Number_of_reads});
}



print R "Coverage\t$coverage\n";
print R "Number_of_Reads_Assembled\t$nreadsAssembled\t$percentReadsAssembled\%\n";
print R "Number_of_Bases_Assembled\t$nbasesAssembled\t$percentBasesAssembled\%\n";

#Print the stats line
print R "$sampleName,";
my $totalInitialReads=($initialStatsR1{Number_of_reads}+$initialStatsR2{Number_of_reads});
print R "$totalInitialReads,";
my $totalStep1Reads=($step1StatsR1{Number_of_reads}+$step1StatsR2{Number_of_reads});
print R "$totalStep1Reads,";
my $percentStep1Reads=getperc($totalStep1Reads,$totalInitialReads);
print R "$percentStep1Reads,";
my $totalStep2Reads=($step2StatsR1{Number_of_reads}+$step2StatsR2{Number_of_reads});
print R "$totalStep2Reads,";
my $percentStep2Reads=getperc($totalStep2Reads,$totalInitialReads);
print R "$percentStep2Reads,";
my $totalStep3Reads=($step3StatsR1{Number_of_reads}+$step3StatsR2{Number_of_reads});
print R "$totalStep3Reads,";
my $percentStep3Reads=getperc($totalStep3Reads,$totalInitialReads);
print R "$percentStep3Reads,";
my $totalStep4Reads=($step4StatsR1{Number_of_reads}+$step4StatsR2{Number_of_reads});
print R "$totalStep4Reads,";
my $percentStep4Reads=getperc($totalStep4Reads,$totalInitialReads);
print R "$percentStep4Reads,";
print R "$nreadsAssembled,";
my $percentReadsAssembled_wrt_initialReads=getperc($nreadsAssembled,$totalInitialReads);
print R "$percentReadsAssembled_wrt_initialReads,";
print R $assemblyStats{"Number_of_contigs"}.",";
print R $assemblyStats{"Largest_contig_size"}.",";
print R $assemblyStats{"Smallest_contig_size"}.",";
print R $assemblyStats{"Average_contig_size"}.",";
print R $assemblyStats{"N50Size"}.",";
print R $assemblyStats{"N50Number"}.",";
print R $assemblyStats{"N50Avg_contig_size"}.",";
print R $assemblyStats{"N90Size"}.",";
print R $assemblyStats{"N90Number"}.",";
print R $assemblyStats{"N90Avg_contig_size"}.",";
print R $assemblyStats{"GenomeSize"}.",";
print R $assemblyStats{"GC"}.",";
print R $assemblyStats{"Number_of_Ns"}.",";
print R $assemblyStats{"Percent_Ns"}.",";
print R $assemblyStats{"Number_of_contigs>10k"}.",";
print R $assemblyStats{"Percent_of_contigs>10k"}.",";
print R $assemblyStats{"Total_size_of_contigs>10k"}.",";
print R $assemblyStats{"Percent_of_genome_in_contigs>10k"}.",";
print R "$coverage\n";


close R;

unless (defined $keepFiles) {
    for my $file (@todelete) {
	unlink $file;
    }
}

exit;


sub usage {
print <<EOF;
MiSeq Assembly Pipeline

Steps:
1. Remove duplicates
2. Quality trimming
3. Low complexity filtering
4. Poly A/T and N filtering
5. Assembly
6. Basic assembly statistics
7. Coverage estimate

Note:
1. Only works with paired end Illumina data.

Author: Vishal N. Koparde, Ph. D.
Created: 140122
Modified: 140122

options:
-1 read1.fastq or read1.fastq.gz
-2 read2.fastq or read2.fastq.gz
-s name of the sample
-n number of cpus
-nodup do not remove duplicates
-nopat do not trim polyA/T from 5' and 3' ends
-nodust do not perform dusting (low complexity filtering)
-noqt do not perform quality trimming
-nocov do not calculate assembly coverage information
-k keep intermediate files
-m minimum read length to consider (default=50)
-w sliding window width for quality trimming (default=9)
-wqt sliding window quality threshold (default=25)
-arqt average read quality threshold (default=30)

EOF
exit 1;
}

sub getReadStats {
    my ($fq)=@_;
    my $fq_stats="${fq}.numreads+";
    my $nreads=`$cfgHash{"fastq_num_reads+"} $fq > $fq_stats`;
    my $tmp=Config::General->new($fq_stats);
    push @todelete,$fq_stats;
    return $tmp->getall;    
}

sub getperc {
    my ($num,$den)=@_;
    my $p=100.0*$num/$den;
    my $q=sprintf "%.2f",$p;
    return $q;
}
