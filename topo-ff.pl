#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: topo-ff.pl
#
#        USAGE: ./topo-ff.pl  TOPOLOGY CHARMM OUTPUT
#
#  DESCRIPTION: Parametrize the TOPOLOGY file for LAMMPS generated with 
#               the VMD plugin 'topotools' with a given CHARMM force field,
#               i.e. use the provided CHARMM topology and parameter file,
#               then write the completed LAMMPS data to a OUTPUT file.
#
#         TODO - add improper
#              - see Axel's notes: NBFIX and water!
#
#        NOTES: Parameters in CHARMM file cannot be separated by empty lines,
#               they have the meaning to separate the different sections.
#               
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

# pragmas
use utf8;
use strict;
use warnings;

# consistency check
my $num_args = $#ARGV + 1;
die("Usage: topo-ff.pl TOPOLOGY CHARMM OUTPUT\n") if $num_args != 3;

# files
my $topo = $ARGV[0];
my $ff  = $ARGV[1];
my $output = $ARGV[2];

open(my $FH_TOPO, '<:encoding(UTF-8)', $topo)
  or die "Could not open file '$topo': $!";

open(my $FH_FF, '<:encoding(UTF-8)', $ff)
  or die "Could not open file '$ff': $!";

open (my $OUT, '>', $output)
    or die "Could not open file '$output': $!";

open (my $NOTFOUND, '>', 'not_found')
    or die "Could not open file 'not_found': $!";

# coefficients
my @masses;
my @pairs;
my @bonds;
my @angles;
my @dihedrals;
my @impropers;

sub handle_pair_coefficients {
    print $NOTFOUND "Pair Coeffs\n\n";
    print $OUT "Pair Coeffs\n\n";
    my $index = 1;
    while (my $line = <$FH_TOPO>) {
        if ($line =~ /^#.*/) {
            $line =~ s/^#\s*\d+\s*(\w+)/$1/;
            my $found = 0;
            for my $pair (@pairs) {
                chomp($line);
                chomp($pair);    
                if ($pair =~ /^$line/) {
                    my @output = ($pair =~ /(\d+\.\d+)/g);
                    my $val1 = abs($output[1]);
                    my $val2 = 2.0*($output[2]*2.0**(-1/6));
                    print $OUT "$index $val1 $val2 $val1 $val2 # $line\n";
                    $index++;
                    $found = 1;
                    last;
                }
            }
            if ($found == 0) {
                print $NOTFOUND $line;
            }
        } else {
            return;
        }
    }
}

sub handle_bond_coefficients {
    print $NOTFOUND "Bond Coeffs\n\n";
    print $OUT "Bond Coeffs\n\n";
    my $index = 1;
    while (my $line = <$FH_TOPO>) {
        if ($line =~ /^#.*/) {
            chomp($line);
            $line =~ s/^#.*?(\w+-\w+)/$1/;
            my @fromto = split /-/, $line;
            my $found = 0;
            for my $bond (@bonds) {
                if ($bond =~ /^$fromto[1]\s+$fromto[0]/) {
                    my @output = ($bond =~ /(\d+\.\d+)/g);
                    print $OUT "$index @output # $line\n";
                    $index++;
                    $found = 1;
                    last; # break OUT of loop if duplicated bond coefficients are specified, take the first definition (all are equivalent)
                }
            }
            if ($found == 0) {
                print $NOTFOUND "Line: $line\n";
            }
        } else {
            return; 
        }
    }
}

sub handle_angle_coefficients {
    print $OUT "Angle Coeffs\n\n";
    print $NOTFOUND "Angle Coeffs\n\n";
    my $index = 1;
    while (my $line = <$FH_TOPO>) {
        if ($line =~ /^#.*/) {
            chomp($line);
            $line =~ s/^#.*?(\w+-\w+-\w+)/$1/;
            my @fromto = split /-/, $line;
            my $found = 0;
            for my $angle (@angles) {
                 if ($angle =~ /$fromto[2]\s*$fromto[1]\s*$fromto[0]/) {
                    my @output = ($angle =~ /(\d+\.\d+)/g);
                    print $OUT "$index @output # $line\n";
                    $found = 1;
                    $index++;
                    last; 
                }
            }
            if ($found == 0) {
                print $NOTFOUND "Line: $line\n";
            }
        } else {
            return; 
        }
    }
}

sub handle_dihedral_coefficients {
    print $NOTFOUND "Dihedral Coeffs\n\n";
    print $OUT "Dihedral Coeffs\n\n";
    my $index = 1;
    while (my $line = <$FH_TOPO>) {
        if ($line =~ /^#.*/) {
            chomp($line);
            $line =~ s/^#.*?(\w+-\w+-\w+-\w+)/$1/;
            my @fromto = split /-/, $line;
            my $found = 0;
            for my $dihedral (@dihedrals) {

                if ($dihedral =~ /$fromto[3]\s*$fromto[2]\s*$fromto[1]\s*$fromto[0]/) {
                    my @output = ($dihedral=~ /(\d+\.\d+)\s*(\d+)\s*(\d+\.\d+)/);
                    print $OUT "$index @output # $line \n";
                    $index++;
                    $found = 1;
                    last; 
                }

                 if ($dihedral =~ /$fromto[0]\s*$fromto[1]\s*$fromto[2]\s*$fromto[3]/) {
                    my @output = ($dihedral=~ /(\d+\.\d+)\s*(\d+)\s*(\d+\.\d+)/);
                    print $OUT "$index @output # $line \n";
                    $index++;
                    $found = 1;
                    last; 
                }
            }

            # try wildcard match
            if ($found == 0) {
                for my $dihedral (@dihedrals) {
                    if ($dihedral =~ /X\s*$fromto[1]\s*$fromto[2]\s*X/) {
                        my @output = ($dihedral=~ /(\d+\.\d+)\s*(\d+)\s*(\d+\.\d+)/);
                        print $OUT "$index @output # $line \n";
                        $index++;
                        $found = 1;
                        last; 
                    }

                    if ($dihedral =~ /X\s*$fromto[2]\s*$fromto[1]\s*X/) {
                        my @output = ($dihedral=~ /(\d+\.\d+)\s*(\d+)\s*(\d+\.\d+)/);
                        print $OUT "$index @output # $line \n";
                        $index++;
                        $found = 1;
                        last; 
                    }
                }

                # no success - need manual fix
                if ($found == 0) {
                    print $NOTFOUND "Line: $line\n";
                }
            }
        } else {
            return; 
        }
    }
}

sub handle_improper_coefficients {
    
}

sub read_coefficients {
    my $type = shift;
    while (<$FH_FF>) {
        chomp;
        if (/^\s*$/) { # empty line marks end of given type section
            last;
        } else {
            if (! (/^!/ || /^\s+!/) ) { # comment line
                push @$type, $_;
            }
        }
    }
}


my %types = ( "BONDS"     => \@bonds,
              "NONBONDED" => \@pairs,
              "ANGLES"    => \@angles,
              "DIHEDRALS" => \@dihedrals,
              "IMPROPERS" => \@impropers
            );

for my $type (keys %types) {
    seek $FH_FF, 0, 0;
    while (my $line = <$FH_FF>) {
        if ($line =~/^$type/) {
            read_coefficients $types{$type};
        }
    }
}

print "Masses: ";
print scalar @masses;
print " ";
print "Bonds: ";
print scalar @bonds;
print " ";
print "Pairs: ";
print scalar @pairs;
print " ";
print "Angles: ";
print scalar @angles;
print " ";
print "Dihedrals: ";
print scalar @dihedrals;
print " ";
print "Impropers: ";
print scalar @impropers;
print "\n";

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
    } elsif ($line =~ /^# Dihedral Coeffs/) { 
         if (my $next_line = <$FH_TOPO>) {
            if ($next_line =~ /^#$/) {
                print $OUT "\n";
                handle_dihedral_coefficients();
                print $OUT "\n";
            }
        }
    } else {
        print $OUT $line;
    }
}

#print join "\n", @dihedrals;
