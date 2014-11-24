#!/usr/bin/env perl

use FindBin;
use lib "$FindBin::RealBin/lib";

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use Schedule::SGELK;

my $dir=`dirname $0`;chomp($dir);
$ENV{PATH}="$ENV{PATH}:$dir/kSNP";

my $sge=Schedule::SGELK->new(-verbose=>1,-numnodes=>5,-numcpus=>8);
my @fastqExt=qw(.fastq .fastq.gz);
my @fastaExt=qw(.fasta .fna .fa);
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
    logmsg "Adding merged fastas to the database";
    addFastaToDb($fasta,$db,$settings);

    logmsg "Adding finished genomes to the database";
    addFinishedGenomesToDb($reads,$db,$settings);
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
  system("cp $db $db.bak"); die if $?;

  logmsg "Now removing one at a time";
  die "Removing multiple genomes at a time is not supported right now" if (@$reads>1);
  for(@$reads){
    my($name,$path,$suffix)=fileparse($_,@fastqExt);
    my $fasta="$name.fasta";

    # find which line the read is on
    my $lineNumber=`grep -n '$fasta' $db.reads | cut -f 1 -d:`; chomp($lineNumber);
    my $fastaLineNumber=$lineNumber * 2 - 1;
    my $line2=$fastaLineNumber+1;
    my $deleteCommand="$fastaLineNumber,$line2"."d";
    my $indexDeleteCommand=$lineNumber."d";

    logmsg "Removing $_ from database";
    $sge->set("jobname","removeFromFasta-$db");
    $sge->pleaseExecute("sed --in-place '$deleteCommand' '$db'");
    for my $index("$db.reads","$db.assemblies"){
      $sge->set("jobname","removeFromIndex-$index");
      $sge->pleaseExecute("sed --in-place '$indexDeleteCommand' $index");
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
  my (@out,@fastaRead);
  for my $r(@$reads){
    my($name,$path,$suffix)=fileparse($r,@fastqExt,@fastaExt);
    my $fasta="$$settings{tempdir}/$name.fasta";
    next if(!$suffix); # meaning, it is not a fastq
    push(@fastaRead,$fasta);
    #die "ERROR: suffix not understood for $r" if(!$suffix);
    logmsg "$r => $fasta";
    next if(-e $fasta);

    my $fastqRegex="(".join("|",@fastqExt).")\$";
    my $fastaRegex="(".join("|",@fastaExt).")\$";
    if($suffix=~/$fastqRegex/){
      # execute the command to convert
      my $command="fastq_to_fasta -Q33";
      if($suffix=~/gz/){
        $command="gunzip -c '$r' | $command";
      } else {
        $command="$command < '$r'";
      }
      # output to a temp file and then move it to the correct name, to ensure that file existance checks work
      $command.=" > $fasta.tmp";
      $command.=" ; mv -v $fasta.tmp $fasta";
      $sge->set("jobname","toFasta-$name");
      $sge->pleaseExecute($command);
    } elsif($suffix=~/$fastaRegex/){
      $sge->pleaseExecute("cp '$r' '$fasta'");
    } else {
      die "ERROR: could not determine the filetype of $r";
    }
  }
  $sge->wrapItUp();

  # merge the fasta reads
  for my $f (@fastaRead){
    my($name,$path,$suffix)=fileparse($f,".fasta");
    my $merged="$$settings{tempdir}/$name.merged.fasta";
    logmsg "$f => $merged";
    push(@out,$merged);
    next if(-e $merged);

    $sge->pleaseExecute("merge_fasta_reads '$f' > $merged.tmp && mv -v $merged.tmp $merged");
  }
  $sge->wrapItUp();
  return \@out;
}

sub addFastaToDb{
  my($fasta,$db,$settings)=@_;
  my $sortedFiles=join(" ",@$fasta);
  my @indexedNames=map{basename($_,".merged.fasta").".fasta"} @$fasta;
  my $indexedNames=join("\n",@indexedNames);
  claimDb($db,$settings);
  $sge->set("jobname","mergingTheDatabase-Reads");
  $sge->pleaseExecute("cat $sortedFiles >> '$db'");
  $sge->set("jobname","addingIndexes");
  $sge->pleaseExecute("echo '$indexedNames' >> '$db.reads'");
  $sge->wrapItUp();
  unlockDb($db,$settings);
}

sub addFinishedGenomesToDb{
  my($fasta,$db,$settings)=@_;

  # merging
  my @merged;
  my @indexedNames;
  for my $f(@$fasta){
    my($name,$path,$suffix)=fileparse($f,@fastaExt);
    my $merged="$$settings{tempdir}/$name.mergedContigs.fasta";
    logmsg "$f => $merged";
    push(@merged,$merged);
    push(@indexedNames,$name);
    next if(-e $merged);
    
    logmsg "merge_fasta_contigs '$f' > $merged.tmp && mv -v $merged.tmp $merged";die;
    die;

    $sge->pleaseExecute("merge_fasta_contigs '$f' > $merged.tmp && mv -v $merged.tmp $merged");
  }
  $sge->wrapItUp();

  # adding to db
  claimDb($db,$settings);
  my $indexedNames=join("\n",@indexedNames);
  $sge->pleaseExecute("cat ".join(" ",@merged)." >> '$db'",{jobname=>"mergingTheDatabase-FinishedGenome"});
  $sge->pleaseExecute("echo '$indexedNames' >> '$db.assemblies",{jobname=>"addingIndexes"});
  $sge->wrapItUp();

  return 1;
}

########################
# database micromanaging
########################
sub claimDb{
  my($db,$settings)=@_;
  while(isLocked($db,$settings)){
    warn "WARNING: Database is locked! Waiting 10 seconds...\n";
    sleep 10;
  }
  #die "ERROR: Database is locked!\n  $db" if(isLocked($db,$settings));
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
  -t tmp/
  "
}
