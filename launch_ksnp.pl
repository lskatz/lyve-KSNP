#!/usr/bin/env perl

use FindBin;
use lib "$FindBin::RealBin/lib";

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use Schedule::SGELK;

my $sge=Schedule::SGELK->new(-verbose=>1,-numnodes=>5,-numcpus=>8);
sub logmsg {local $0=basename $0;$|++;my $FH = *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n"; $|--;}
exit(main());
sub main{
  my $settings={};
  GetOptions($settings,qw(ref=s@ help tempdir=s));

  my $ref=$$settings{ref} || die usage();
  my @reads=@ARGV;
  die usage() if (@reads<1);
  $$settings{tempdir}||=die usage();
  $$settings{numcpus}||=8;
  $$settings{numnodes}||=5;
  $sge->set("workingdir",$$settings{tempdir});
  for (qw(workingdir numnodes numcpus)){
    $sge->set($_,$$settings{$_});
  }
  
  # check if ksnp is in the path
  system("which kSNP >& /dev/null");
  die "ERROR: could not find kSNP in the path" if $?;
  
  logmsg "Merging assemblies into one file";
  my $asmMerged=mergeAssemblies($ref,$settings);
  logmsg "Converting fastq reads to fasta";
  my @fasta=fastqToFasta(\@reads,$settings);
  logmsg "Merging raw fasta reads";
  my $allMerged=mergeReads_withAsm(\@fasta,$asmMerged,$settings);
  logmsg "Running ksnp on $allMerged";
  my $outDir=ksnp($allMerged,$settings);
  logmsg "All files are in $outDir";

  return 0;
}

sub mergeAssemblies{
  my($ref,$settings)=@_;
  my $mergedFile="$$settings{tempdir}/ref.merged.fasta";
  return $mergedFile if(-e $mergedFile && -s $mergedFile > 1);
  $sge->pleaseExecute_andWait("cat ".join(" ",@$ref)." > $$settings{tempdir}/ref.fasta",$settings);
  $sge->pleaseExecute_andWait("merge_fasta_contigs $$settings{tempdir}/ref.fasta > $mergedFile",$settings);
  die "ERROR: merged file $mergedFile does not contain an assembly!" if(!-e $mergedFile || -s $mergedFile < 1);
  return $mergedFile;
}

sub fastqToFasta{
  my($reads,$settings)=@_;
  my @out;
  print "\n";
  for my $r(@$reads){
    my $b=basename($r,qw(.fastq));
    my $out="$$settings{tempdir}/$b.fasta";
    push(@out,$out);
    #logmsg "$r => $out";
    next if (-e $out);
    open(FASTQ,$r) or die "ERROR: could not open $r for reading:$!";
    open(OUT,">$out.tmp") or die "ERROR: could not open $out for writing: $!";
    my $i=0;
    while(<FASTQ>){
      $i++;
      my $mod=$i%4;
      if($mod eq 1){
        s/^@/>/;
        print OUT $_;
      } elsif($mod eq 2){
        print OUT $_;
      }
    }
    close OUT;
    close FASTQ;
    system("mv $out.tmp $out"); die if $?;
    print "."; # add a dot to show progress
  }
  print "\n";
  return @out;
}

sub mergeReads_withAsm{
  my($fasta,$mergedAsm,$settings)=@_;
  my(@merged,@cmd);
  my $mergedOut="$$settings{tempdir}/reads.merged.fasta";
  for(@$fasta){
    my $b=basename($_,qw(.fasta));
    my $merged="$$settings{tempdir}/$b.merged.fasta";
    #logmsg "Merging $_ => $merged";
    push(@cmd,"merge_fasta_reads $_ > $merged") if(!-e $merged || -s $merged <1);
    push(@merged,$merged);
  }
  $sge->pleaseExecute_andWait(\@cmd,$settings);

  logmsg "Concatenating $mergedAsm with all merged raw read files => $mergedOut";
  $sge->pleaseExecute_andWait([
    "cat $mergedAsm ".join(" ",@merged)." > $mergedOut",
    "genome_names $mergedAsm > $$settings{tempdir}/genomes.finished",
  ],$settings) if(!-e $mergedOut || -s $mergedOut <1);
  open(UNFINISHED,">$$settings{tempdir}/genomes.unfinished") or die "Could not open $$settings{tempdir}/genomes.unfinished: $!";
  for(@$fasta){
    print UNFINISHED basename($_)."\n";
  }
  close UNFINISHED;
  return $mergedOut;
}
  
sub ksnp{
  my($mergedFasta,$settings)=@_;
  my $finishedList="$$settings{tempdir}/genomes.finished";
  my $unfinishedList="$$settings{tempdir}/genomes.unfinished";
  my $outDir="ksnpOut";
  $sge->pleaseExecute_andWait("kSNP -f $mergedFasta -k 25 -d $outDir -p $finishedList -u $unfinishedList",$settings);
  return $outDir;
}

sub usage{
  "Runs kSNP from a set of cleaned fastq files
  Usage: $0 *.fastq -ref reference.fasta [-ref reference2.fasta] -t tempdir/
  "
}
