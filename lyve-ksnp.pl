#!/usr/bin/env perl

use FindBin;
use lib "$FindBin::RealBin/lib";

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use Schedule::SGELK;
use Filelock qw/:all/;

my $dir=`dirname $0`;chomp($dir);
$ENV{PATH}="$ENV{PATH}:$dir/kSNP";

my $sge;
sub logmsg {local $0=basename $0;$|++;my $FH = *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n"; $|--;}
exit(main());
sub main{
  my $settings={};
  GetOptions($settings,qw(db=s help tempdir=s outdir=s kmerlength=i keep));

  my $ref=$$settings{ref} || "";
  $$settings{tempdir}||=die "ERROR: need temporary directory!\n". usage();
  $$settings{db}||=die "ERROR: need fasta database!\n".usage();
  $$settings{numcpus}||=8;
  $$settings{numnodes}||=5;
  $$settings{outdir}||="ksnpOut";
  $$settings{kmerlength}||=25;

  $sge=Schedule::SGELK->new(-verbose=>1,-numnodes=>5,-numcpus=>8);
  $sge->set('qsubxopts',$$settings{qsubxopts}) if($$settings{qsubxopts});
  $sge->set("workingdir",$$settings{tempdir});
  for (qw(workingdir numnodes numcpus)){
    $sge->set($_,$$settings{$_});
  }
  
  # check if ksnp is in the path
  system("which kSNP >& /dev/null");
  die "ERROR: could not find kSNP in the path" if $?;
  
  logmsg "Running ksnp on $$settings{db}";
  my $outDir=ksnp($$settings{db},$settings);
  logmsg "All files are in $outDir";

  return 0;
}

sub mergeAssemblies{
  my($ref,$settings)=@_;
  my $mergedFile="$$settings{tempdir}/ref.merged.fasta";
  system("touch $mergedFile") if(!$ref);
  return $mergedFile if(-e $mergedFile);
  $sge->pleaseExecute_andWait("cat ".join(" ",@$ref)." > $$settings{tempdir}/ref.fasta",$settings);
  $sge->pleaseExecute_andWait("merge_fasta_contigs $$settings{tempdir}/ref.fasta > $mergedFile",$settings);
  die "ERROR: merged file $mergedFile does not contain an assembly!" if(!-e $mergedFile || -s $mergedFile < 1);
  return $mergedFile;
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
  my $outDir=$$settings{outdir};
  $sge->pleaseExecute_andWait("kSNP_notrees -f $mergedFasta -k $$settings{kmerlength} -d $outDir -p $finishedList -u $unfinishedList -n $$settings{numcpus}",$settings);
  # remove temporary directory in $outdir called TemporaryFilesToDelete
  system("rm -rfv $outDir/TemporaryFilesToDelete") unless($$settings{keep});
  return $outDir;
}

sub usage{
  "Runs kSNP from a 'kSNP database'. You can make this db from lyve-manage-ksnp.pl
  Usage: $0 -d database.fasta -o outdir/ -t tmpdir/
  -kmerlength kmer length (default: 25)
  -keep to keep temporary files
  "
}
