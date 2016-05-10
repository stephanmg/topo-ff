#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: topo-ff.pl
#
#        USAGE: ./topo-ff.pl  TOPOLOGY CHARMM OUTPUT
#
#  DESCRIPTION: Parametrize the TOPOLOGY file for LAMMPS generated with 
#               the VMD plugin 'topotools' with a given CHARMM force field 
#               and write the completed LAMMPS data to a OUTPUT file.
#
#         TODO add improper, angles and dihedrals coefficients/see notes Axel!
#              ignore lines startign with ! indicating a comment
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: None so far.
#        NOTES: In case of any questions send me an email.
#       AUTHOR: Stephan Grein (grein@temple.edu)
# ORGANIZATION: Temple University
#      VERSION: 1.0
#      CREATED: 05/03/16 15:51:36
#     REVISION: cf. git log 
#===============================================================================

use utf8;
use strict;
use warnings;

my $num_args = $#ARGV + 1;
die("Usage: topo-ff.pl TOPOLOGY CHARMM OUTPUT\n") if $num_args != 3;

my $topo = $ARGV[0];
my $ff  = $ARGV[1];
my $output = $ARGV[2];

open(my $FH_TOPO, '<:encoding(UTF-8)', $topo)
  or die "Could not open file '$topo': $!";

open(my $FH_FF, '<:encoding(UTF-8)', $ff)
  or die "Could not open file '$ff': $!";

open (my $OUT, '>', $output)
    or die "Could not open file '$output': $!";

my @pairs;
my @bonds;
my @angles;
my @dihedrals;
my @impropers;

=begin
while (my $line = <$FH_FF>) {
    # TODO 5 if statements can be compressed into a for loop over [pairs, bonds, angles, dihedrals, impropers]
    if ($line =~ /^Pair Coeffs/) {
        if (my $next_line = <$FH_FF>) {
            if ($next_line =~ /^$/) {
                while (my $pair = <$FH_FF>) {
                    if ($pair =~ /\d+.*#.*/) {
                        push @pairs, $pair;
                        } else {
                        last;
                    }
                }
            }
        }
    }
 if ($line =~ /^Bond Coeffs/) {
        if (my $next_line = <$FH_FF>) {
            if ($next_line =~ /^$/) {
                while (my $bond = <$FH_FF>) {
                    if ($bond =~ /\d+.*#.*/) {
                        push @bonds, $bond;
                        } else {
                        last;
                    }
                }
            }
        }
    }
  if ($line =~ /^Angle Coeffs/) {
        if (my $next_line = <$FH_FF>) {
            if ($next_line =~ /^$/) {
                while (my $angle = <$FH_FF>) {
                    if ($angle =~ /\d+.*#.*/) {
                        push @angles, $angle;
                        } else {
                        last;
                    }
                }
            }
        }
    }
  if ($line =~ /^Dihedral Coeffs/) {
        if (my $next_line = <$FH_FF>) {
            if ($next_line =~ /^$/) {
                while (my $dihedral = <$FH_FF>) {
                    if ($dihedral =~ /\d+.*#.*/) {
                        push @dihedrals, $dihedral;
                        } else {
                        last;
                    }
                }
            }
        }
    }
 if ($line =~ /^Improper Coeffs/) {
        if (my $next_line = <$FH_FF>) {
            if ($next_line =~ /^$/) {
                while (my $improper = <$FH_FF>) {
                    if ($improper =~ /\d+.*#.*/) {
                        push @impropers, $improper;
                        } else {
                        last;
                    }
                }
            }
        }
    }
}

sub handle_angle_coeffs {
    print $OUT "Angle Coeffs\n\n";
    my $index = 1;
    while (my $line = <$FH_TOPO>) {
        if ($line =~ /^#.*/) {
            chomp($line);
            $line =~ s/^#.*?(\w+-\w+-\w+)/$1/;
            my @fromto = split /-/, $line;
            for my $angle (@angles) {
                if ($angle =~ /$fromto[0]\s*$fromto[1]\s*$fromto[2]/) {
                    $angle =~ s/\d+\s*(.*)#.*/$1/;
                    print $OUT "$index $angle";
                    $index++;
                    last; 
                }
 if ($angle =~ /$fromto[1]\s*$fromto[0]\s*$fromto[2]/) {
                    $angle =~ s/\d+\s*(.*)#.*/$1/;
                    print $OUT "$index $angle";
                    $index++;
                    last; 
                }

            }
        } else {
            return; 
        }
    }
}

sub handle_dihedral_coeffs {
    print $OUT "Dihedral Coeffs\n\n";
    my $index = 1;
    while (my $line = <$FH_TOPO>) {
        if ($line =~ /^#.*/) {
            chomp($line);
            $line =~ s/^#.*?(\w+-\w+-\w+-\w+)/$1/;
            my @fromto = split /-/, $line;
            for my $dihedral (@dihedrals) {
                if ($dihedral =~ /$fromto[0]\s*$fromto[1]\s*$fromto[2]\s*$fromto[3]/) {
                    my $OUTput = $dihedral;
                    $OUTput =~ s/\d+\s*(.*)#.*/$1/;
                    print $OUT "$index $OUTput";
                    $index++;
                    last; 
                }

               if ($dihedral =~ /$fromto[1]\s*$fromto[0]\s*$fromto[3]\s*$fromto[2]/) {
                    my $OUTput = $dihedral;
                    $OUTput  =~ s/\d+\s*(.*)#.*/$1/;
                    print $OUT "$index $OUTput";
                    $index++;
                    last; 
                }

                if ($dihedral =~ /$fromto[3]\s*$fromto[2]\s*$fromto[1]\s*$fromto[0]/) {
                    my $OUTput = $dihedral;
                    $OUTput  =~ s/\d+\s*(.*)#.*/$1/;
                    print $OUT "$index $OUTput";
                    $index++;
                    last; 
                }

                if ($dihedral =~ /$fromto[2]\s*$fromto[1]\s*$fromto[3]\s*$fromto[0]/) {
                    my $OUTput = $dihedral;
                    $OUTput =~ s/\d+\s*(.*)#.*/$1/;
                    print $OUT "$index $OUTput";
                    $index++;
                    last; 
                }
            }
        } else {
            return; 
        }
    }
}

sub handle_improper_coeffs {
    print $OUT "Improper Coeffs\n\n";
    my $index = 1;
    while (my $line = <$FH_TOPO>) {
        if ($line =~ /^#.*/) {
            chomp($line);
            $line =~ s/^#.*?(\w+-\w+-\w+-\w+)/$1/;
            my @fromto = split /-/, $line;
            for my $improper (@impropers) {
                if ($improper =~ /$fromto[2]\s*$fromto[1]\s*$fromto[3]\s*$fromto[0]/) {
                    my $OUTput = $improper;
                    $improper =~ s/\d+\s*(.*)#.*/$1/;
                    print $OUT "$index $improper";
                    $index++;
                    last; 
                }

                 if ($improper =~ /$fromto[0]\s*$fromto[1]\s*$fromto[2]\s*$fromto[3]/) {
                    my $OUTput = $improper;
                    $improper =~ s/\d+\s*(.*)#.*/$1/;
                    print $OUT "$index $improper";
                    $index++;
                    last; 
                }
    if ($improper =~ /$fromto[3]\s*$fromto[1]\s*$fromto[0]\s*$fromto[2]/) {
                    my $OUTput = $improper;
                    $improper =~ s/\d+\s*(.*)#.*/$1/;
                    print $OUT "$index $improper";
                    $index++;
                    last; 
                }



            }
        } else {
            return; 
        }
    }
}
=end
=cut

sub handle_pair_coefficients {
    print $OUT "Pair Coeffs\n\n";
    my $index = 1;
    while (my $line = <$FH_TOPO>) {
        if ($line =~ /^#.*/) {
            $line =~ s/^#\s*\d+\s*(\w+)/$1/;
            for my $pair (@pairs) {
                chomp($line);
                chomp($pair);    
                if ($pair =~ /^$line/) {
                    my @output = ($pair =~ /(\d+\.\d+)/g);
                    my $val1 = abs($output[1]);
                    my $val2 = 2.0*($output[2]*2.0**(-1/6));
                    print $OUT "$index $val1 $val2 $val1 $val2 # $line\n";
                    $index++;
                    last;
                }
            }
        } else {
            return;
        }
    }
}

sub handle_bond_coefficients {
    print $OUT "Bond Coeffs\n\n";
    my $index = 1;
    while (my $line = <$FH_TOPO>) {
        if ($line =~ /^#.*/) {
            chomp($line);
            $line =~ s/^#.*?(\w+-\w+)/$1/;
            my @fromto = split /-/, $line;
            for my $bond (@bonds) {
                if ($bond =~ /^$fromto[1]\s+$fromto[0]/) {
                    my @output = ($bond =~ /(\d+\.\d+)/g);
                    print $OUT "$index @output # $line\n";
                    $index++;
                    last; # break OUT of loop if duplicated bond coefficients are specified, take the first definition (all are equivalent)
                }
            }
            
        } else {
            return; 
        }
    }
}

sub handle_angle_coefficients {
    print $OUT "Angle Coeffs\n\n";
    my $index = 1;
    while (my $line = <$FH_TOPO>) {
        if ($line =~ /^#.*/) {
            chomp($line);
            $line =~ s/^#.*?(\w+-\w+-\w+)/$1/;
            my @fromto = split /-/, $line;
            for my $angle (@angles) {
                 if ($angle =~ /$fromto[2]\s*$fromto[1]\s*$fromto[0]/) {
                    my @output = ($angle =~ /(\d+\.\d+)/g);
                    print $OUT "$index @output # $line\n";
                    $index++;
                    last; 
                }

            }
        } else {
            return; 
        }
    }
}




sub read_coefficients {
    my $type = shift;
    while (<$FH_FF>) {
        chomp;
        if (/^\s*$/) { # empty line marks end of given type section
            last;
        } else {
            if (! (/^!/ || /^\s+!/) ) {
                push @$type, $_;
            }
        }
    }
}

my %types = ( "BONDS"     => \@bonds,
              "NONBONDED" => \@pairs,
              "ANGLES"    => \@angles
            );

for my $type (keys %types) {
    seek $FH_FF, 0, 0;
    while (my $line = <$FH_FF>) {
        if ($line =~/^$type/) {
            read_coefficients $types{$type};
        }
    }
}

print "Bonds: ";
print scalar @bonds;
print " ";
print "Pairs: ";
print scalar @pairs;
print " ";
print "Angles: ";
print scalar @angles;

while (my $line = <$FH_TOPO>) {
    if ($line =~ /^# Pair Coeffs/) {
        if (my $next_line = <$FH_TOPO>) {
            if ($next_line =~ /^#$/) {
                print $OUT "\n";
                handle_pair_coefficients();
                print $OUT "\n";
            }
        }
    } elsif ($line =~ /^# Bond Coeffs/) {
        if (my $next_line = <$FH_TOPO>) {
            if ($next_line =~ /^#$/) {
                print $OUT "\n";
                handle_bond_coefficients();
                print $OUT "\n";
            }
        }
    } elsif ($line =~ /^# Angle Coeffs/) { 
        if (my $next_line = <$FH_TOPO>) {
            if ($next_line =~ /^#$/) {
                print $OUT "\n";
                handle_angle_coefficients();
                print $OUT "\n";
            }
        }
    } else {
        print $OUT $line;
    }
}
