#!/usr/bin/env perl

use strict;
use warnings;

# simultaneously lock the database and check if it is already locked.
sub claimDb{
  my($db,$settings)=@_;
  die "ERROR: Database is locked!\n  $db" if(isLocked($db,$settings));
  lockDb($db,$settings);
  return 1;
}
# unlock the database
sub unlockDb{
  my($db,$settings)=@_;
  unlink lockFilename($db,$settings);
  return 1;
}
# check if the database is locked
sub isLocked{
  my($db,$settings)=@_;
  my $lockfilename=lockFilename($db,$settings);
  return 1 if(-e $lockfilename);
  return 0;
}

##################
# helper functions
##################
# lock the database
sub lockDb{
  my($db,$settings)=@_;
  my $lockfilename=lockFilename($db,$settings);
  system("touch '$lockfilename'"); die if $?;
  return 1;
}
# create the database lockfile name
sub lockFilename{
  my($db,$settings)=@_;
  return ".$db.locked";
}

1;
