#!/usr/bin/env perl

use strict;
use warnings;

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

1;
