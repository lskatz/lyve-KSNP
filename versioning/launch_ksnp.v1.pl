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
  pleaseExecute_andWait("cat ".join(" ",@$ref)." > $$settings{tempdir}/ref.fasta",$settings);
  pleaseExecute_andWait("merge_fasta_contigs $$settings{tempdir}/ref.fasta > $mergedFile",$settings);
  die "ERROR: merged file $mergedFile does not contain an assembly!" if(!-e $mergedFile || -s $mergedFile < 1);
  return $mergedFile;
}

sub fastqToFasta{
  my($reads,$settings)=@_;
  my @out;
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
  }
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
  pleaseExecute_andWait(\@cmd,$settings);

  logmsg "Concatenating $mergedAsm with all merged raw read files => $mergedOut";
  pleaseExecute_andWait([
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
  pleaseExecute_andWait("kSNP -f $mergedFasta -k 25 -d $outDir -p $finishedList -u $unfinishedList",$settings);
  return $outDir;
}

sub pleaseExecute_andWait{
  my($cmd,$settings)=@_;
  my %settings=%$settings; # don't overwrite settings by reference
  $cmd=[$cmd] if(ref($cmd) eq ""); # let cmd be a string but turn it into a list internally
  $settings{mustfinish}=0;
  my(@jobid);
  for(@$cmd){
    my $jobid=pleaseExecute($_,\%settings);
    push(@jobid,$jobid);
    waitOnJobs(\@jobid,\%settings);
  }
  $settings{mustfinish}=1;
  waitOnJobs(\@jobid,\%settings);
}

# submit to the cluster
sub pleaseExecute{
  my($cmd,$settings)=@_;
  local $0=basename $0;
  my $jobid=-1; # -1 is an error state
  $$settings{jobname}||="run$0";
  $$settings{logfile}||="$0.log";
  $$settings{numcpus}||=1;
  return 0 if($cmd=~/^\s*$/); # if there's no command, then no worries

  # Removing the redirect in the command before putting it into qsub (very rudimentary)
  $cmd=~s/ 1?>(\S+)/ /;
  my $outfile=$1 || $$settings{logfile};

  $cmd="echo '$cmd' | qsub -pe smp $$settings{numcpus} -N $$settings{jobname} -cwd -V -o $outfile -e $$settings{logfile} 2>&1";
  my $out=`$cmd`;
  print $out;
  if($out=~/Your job (\d+)/){
    $jobid=$1;
  } else {
    logmsg "WARNING: the last job submitted did not have an obvious jobid. It can't be tracked!";
  }
  return $jobid;
}

# return the job status of the id
sub checkJob{
  my($jobid,$settings)=@_;
  my $state=0;
  open(QSTAT,"qstat|") or die "ERROR: could not execute qstat!";
  while(<QSTAT>){
    s/^\s+|\s+$//g;
    my @F=split /\s+/;
    if($F[0] eq $jobid){
      $state=$F[4];
    }
  }
  close QSTAT;
  return $state;
}

sub usage{
  "Runs kSNP from a set of cleaned fastq files
  Usage: $0 *.fastq -ref reference.fasta [-ref reference2.fasta] -t tempdir/
  "
}

sub waitOnJobs{
  my($jobid,$settings)=@_;
  logmsg "We have reached node capacity ($$settings{numnodes})! Waiting for a job to finish." if(@$jobid >= $$settings{numnodes});
  while(@$jobid > 0){
    for(my $i=0;$i<@$jobid;$i++){
      my $state=checkJob($$jobid[$i],$settings);
      if($state eq 0){
        logmsg "A job finished: $$jobid[$i]";
        splice(@$jobid,$i,1);
        last;
      }
    }
    sleep 1;
    # break out if you don't have to finish yet but you can still add in another job
    last if(!$$settings{mustfinish} && @$jobid<$$settings{numnodes});
  }
  return @$jobid;
}



