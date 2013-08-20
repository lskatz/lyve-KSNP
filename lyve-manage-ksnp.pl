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
my @fastqExt=qw(.fastq .fastq.gz);
sub logmsg {local $0=basename $0;$|++;my $FH = *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n"; $|--;}
exit(main());
sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s db=s action=s));

  my @reads=@ARGV;
  die usage() if (@reads<1);
  $$settings{tempdir}||=die "ERROR: need temporary directory\n". usage();
  $$settings{numcpus}||=1;
  $$settings{numnodes}||=20;
  $$settings{db}||=die "ERROR: need database fasta file\n".usage();
  $$settings{action}||="add";
  $sge->set("workingdir",$$settings{tempdir});
  for (qw(workingdir numnodes numcpus)){
    $sge->set($_,$$settings{$_});
  }
  
  # check if ksnp is in the path
  system("which kSNP >& /dev/null");
  die "ERROR: could not find kSNP in the path" if $?;

  my $db=$$settings{db};
  my $action=lc($$settings{action});
  if($action eq "add"){
    my $reads=ignoreWhatWeAlreadyHave(\@reads,["$db.reads","$db.assemblies"],$settings);
    logmsg "Converting fastq reads to fasta";
    my $fasta=fastqToFasta($reads,$settings);
    logmsg "Adding fasta to database";
    addFastaToDb($fasta,$db,$settings);
  }elsif($action eq "remove"){
    removeSequences(\@reads,$db,$settings);
  }else{
    logmsg "I'm unsure what to do with the action $action !".usage();
    die;
  }

  return 0;
}

sub removeSequences{
  my($reads,$db,$settings)=@_;
  claimDb($db,$settings);
  # back up the database one time instead of many times with sed
  logmsg "Backing up the database to $db.bak";
  #system("cp $db $db.bak"); die if $?;
  for(@$reads){
    my($name,$path,$suffix)=fileparse($_,@fastqExt);
    my $fasta="$name.fasta";
    logmsg "Removing $_ from database";
    $sge->pleaseExecute("sed --in-place '/$fasta/d' $db");
    for my $index("$db.reads","$db.assemblies"){
      $sge->pleaseExecute("sed --in-place '/$fasta/d' $index");
    }
    $sge->wrapItUp();
  }
  unlockDb($db,$settings);
  return 1;
}


sub ignoreWhatWeAlreadyHave{
  my($reads,$index,$settings)=@_;
  my %alreadyHave=();
  system("touch ".join(" ",@$index)); die if $?;
  for (@$index){
    open(INDEX,$_) or die "Could not read file $_:$!";
    while(<INDEX>){
      chomp;
      $alreadyHave{$_}=1;
    }
    close INDEX;
  }
  for(@$reads){
    my($name,$path,$suffix)=fileparse($_,@fastqExt);
    if($alreadyHave{"$name.fasta"}){
      $_="" if($alreadyHave{"$name.fasta"});
      logmsg "Warning: $name is already present in the database. Ignoring.";
    }
  }
  $reads=[grep(!/^\s*$/,@$reads)];
  return $reads;
}

sub fastqToFasta{
  my($reads,$settings)=@_;
  my @out;
  for my $r(@$reads){
    my($name,$path,$suffix)=fileparse($r,@fastqExt);
    my $out="$$settings{tempdir}/$name.fasta";
    push(@out,$out);
    die "ERROR: suffix not understood for $r" if(!$suffix);
    logmsg "$r => $out";

    next if(-e $out);
    # execute the command to convert
    my $command="fastq_to_fasta -Q33";
    if($suffix=~/gz/){
      $command="gunzip -c '$r' | $command";
    } else {
      $command="$command < '$r'";
    }
    # output to a temp file and then move it to the correct name, to ensure that file existance checks work
    $command.=" > $out.tmp";
    $command.=" ; mv -v $out.tmp $out";
    $sge->set("jobname","toFasta-$name");
    $sge->pleaseExecute($command);
  }
  $sge->wrapItUp();
  return \@out;
}

sub addFastaToDb{
  my($fasta,$db,$settings)=@_;
  claimDb($db,$settings);
  # add one fasta at a time, so that the database doesn't get corrupted
  for my $f (@$fasta){
    my($name,$path,$suffix)=fileparse($f);
    $sge->set("jobname","merging-$name");
    $sge->pleaseExecute_andWait("merge_fasta_reads '$f' >> '$db'");
    $sge->pleaseExecute_andWait("echo '$name' >> '$db.reads'");
  }
  $sge->wrapItUp();
  unlockDb($db,$settings);
}

########################
# database micromanaging
########################
sub claimDb{
  my($db,$settings)=@_;
  die "ERROR: Database is locked!\n  $db" if(isLocked($db,$settings));
  lockDb($db,$settings);
  return 1;
}
sub lockDb{
  my($db,$settings)=@_;
  system("touch '$db.locked'"); die if $?;
  return 1;
}
sub unlockDb{
  my($db,$settings)=@_;
  unlink "$db.locked";
  return 1;
}
sub isLocked{
  my($db,$settings)=@_;
  return 1 if(-e "$db.locked");
  return 0;
}

sub usage{
  "Adds raw reads to a database for kSNP
  Usage: $0 *.fastq[.gz] -d database.fasta
  -d database of merged fastas, produced by ksnp executable merge_fasta_reads
  --action indicates one of the following: add, remove. Default: add
  "
}
