#!/usr/bin/env perl 
use 5.000;
use utf8;
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

# cli options
my $topo;
my $ff; 
my $output;
my $verbose = 0;
my $not_found = "not_found.txt";

GetOptions(
    'topology=s' => \$topo,
    'ff=s'       => \$ff,
    'output=s'   => \$output,
    'verbose'  => \$verbose,
    'usage'    => sub { pod2usage(2) },
    'help'     => sub { pod2usage(1) },
    'man'      => sub { pod2usage(-existatus => 0, -verbose => 2) },
);

pod2usage(2) unless defined($topo) and defined($ff) and defined($output);

# open files
open(my $FH_TOPO, "<:encoding(UTF-8)", "$topo")
  or die "Could not open toplogy file '$topo': $!";
open(my $FH_FF, "<:encoding(UTF-8)", "$ff")
  or die "Could not open force field file '$ff': $!";
open (my $OUT, ">:encoding(UTF-8)", "${output}.data")
    or die "Could not open output file '$output': $!";
open (my $CONFIG, ">:encoding(UTF-8)", "${output}.in")
    or die "Could not open config file '${output}.in': $!"; 
open (my $NOTFOUND, ">:encoding(UTF-8)", "$not_found")
    or die "Could not open file '$not_found': $!";

# coefficients
my @pairs;
my @bonds;
my @angles;
my @dihedrals;
my @impropers;

# handle pair coefficients
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

# handle bond coefficients
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
                    last; # break OUT of loop if duplicated bond coefficients are specified, take the first definition 
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

# handle angle coefficients
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
                    if ($#output+1 == 2) {
                       print $OUT "$index @output 0.0 0.0 # $line\n";
                    } else { # Urey-Bradley term
                       print $OUT "$index @output # $line\n";
                    }
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

# handle dihedral coefficients
sub handle_dihedral_coefficients {
    print $NOTFOUND "Dihedral Coeffs\n\n";
    print $OUT "Dihedral Coeffs\n\n";
    my $weighting = 0; # for now no LJ 1-4 weighting, see TODO
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
                    print $OUT "$index @output $weighting # $line \n";
                    $index++;
                    $found = 1;
                    last; 
                }

                 if ($dihedral =~ /$fromto[0]\s*$fromto[1]\s*$fromto[2]\s*$fromto[3]/) {
                    my @output = ($dihedral=~ /(\d+\.\d+)\s*(\d+)\s*(\d+\.\d+)/);
                    print $OUT "$index @output $weighting # $line \n";
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
                        print $OUT "$index @output $weighting # $line \n";
                        $index++;
                        $found = 1;
                        last; 
                    }

                    if ($dihedral =~ /X\s*$fromto[2]\s*$fromto[1]\s*X/) {
                        my @output = ($dihedral=~ /(\d+\.\d+)\s*(\d+)\s*(\d+\.\d+)/);
                        print $OUT "$index @output $weighting # $line \n";
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

# handle improper coefficients
sub handle_improper_coefficients {
    print $OUT "Improper Coeffs\n\n";
    print $NOTFOUND "Improper Coeffs\n\n";
    my $index = 1;
    while (my $line = <$FH_TOPO>) {
        if ($line =~ /^#.*/) {
            chomp($line);
            $line =~ s/^#.*?(\w+-\w+-\w+-\w+)/$1/;
            my @fromto = split /-/, $line;
            my $found = 0;
            for my $improper (@impropers) {
                if ($improper =~ /$fromto[0]\s*$fromto[1]\s*$fromto[2]\s*$fromto[3]/) {
                     my @output = ($improper =~ /(\d+\.\d+)\s*(\d+)\s*(\d+\.\d+)/);
                     # second column is ignored, see CHARMM parameter file
                     print $OUT "$index $output[0] $output[2] # $line \n"; 
                     $index++;
                     $found = 1;
                     last; 
                }

                if ($improper =~ /$fromto[3]\s*$fromto[2]\s*$fromto[1]\s*$fromto[0]/) {
                     my @output = ($improper =~ /(\d+\.\d+)\s*(\d+)\s*(\d+\.\d+)/);
                     print $OUT "$index $output[0] $output[2] # $line \n";
                     $index++;
                     $found = 1;
                     last; 
                }
            }

            # try wildcard match
            if ($found == 0) {
               for my $improper (@impropers) {
                  if ($improper =~ /X\s*$fromto[1]\s*$fromto[2]\s*X/) {
                     my @output = ($improper=~ /(\d+\.\d+)\s*(\d+)\s*(\d+\.\d+)/);
                     print $OUT "$index $output[0] $output[2] # $line \n";
                     $index++;
                     $found = 1;
                     last; 
                  }

                 if ($improper =~ /X\s*$fromto[2]\s*$fromto[1]\s*X/) {
                     my @output = ($improper=~ /(\d+\.\d+)\s*(\d+)\s*(\d+\.\d+)/);
                     print $OUT "$index $output[0] $output[2] # $line \n";
                     $index++;
                     $found = 1;
                     last; 
                 }

                 if ($improper =~ /$fromto[3]\s*X\s*X\s*$fromto[0]/) {
                     my @output = ($improper =~ /(\d+\.\d+)\s*(\d+)\s*(\d+\.\d+)/);
                     print $OUT "$index $output[0] $output[2] # $line \n";
                     $index++;
                     $found = 1;
                     last; 
                }

                 if ($improper =~ /$fromto[0]\s*X\s*X\s*$fromto[3]/) {
                     my @output = ($improper =~ /(\d+\.\d+)\s*(\d+)\s*(\d+\.\d+)/);
                     print $OUT "$index @output # $line \n";
                     $index++;
                     $found = 1;
                     last; 
                }
           }
        }
             # no success - need manual fix
             if ($found == 0) {
                 print $NOTFOUND "Line: $line\n";
             }
        }
    }
}

# helper method to read coefficients for each type
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

# sections in CHARMM param file
my %types = ( "BONDS"     => \@bonds,
              "NONBONDED" => \@pairs,
              "ANGLES"    => \@angles,
              "DIHEDRALS" => \@dihedrals,
              "IMPROPER" => \@impropers
            );

# read in all coefficients from CHARMM param file
for my $type (keys %types) {
    seek $FH_FF, 0, 0;
    while (my $line = <$FH_FF>) {
        if ($line =~/^$type/) {
            read_coefficients $types{$type};
        }
    }
}

# internal statistics (coefficients from provided CHARMM param file)
if ($verbose) {
    print "****************************\n";
    print "*** \tBonds: ";
    print scalar @bonds . "\t ***\n";
    print "*** \tPairs: ";
    print scalar @pairs . "\t *** \n";
    print "*** \tAngles: ";
    print scalar @angles . "\t *** \n";
    print "*** \tDihedrals: ";
    print scalar @dihedrals . "\t *** \n";
    print "*** \tImpropers: ";
    print scalar @impropers . "\t *** \n";
    print "****************************\n";
}

# write data file
while (my $line = <$FH_TOPO>) {
    my $next_line;
    if ($line =~ /^# Pair Coeffs/) {
        if (defined ($next_line = <$FH_TOPO>)) {
           if ($next_line =~ /^#$/) {
                print $OUT "\n";
                handle_pair_coefficients();
                print $OUT "\n";
           }
        }
    } elsif ($line =~ /^# Bond Coeffs/) {
        if (defined($next_line = <$FH_TOPO>)) {
            if ($next_line =~ /^#$/) {
                print $OUT "\n";
                handle_bond_coefficients();
                print $OUT "\n";
            }
        }
    } elsif ($line =~ /^# Angle Coeffs/) { 
        if (defined($next_line = <$FH_TOPO>)) {
            if ($next_line =~ /^#$/) {
                print $OUT "\n";
                handle_angle_coefficients();
                print $OUT "\n";
            }
        }
    } elsif ($line =~ /^# Dihedral Coeffs/) { 
        if (defined($next_line = <$FH_TOPO>)) {
            if ($next_line =~ /^#$/) {
                print $OUT "\n";
                handle_dihedral_coefficients();
                print $OUT "\n";
            }
        }
    } elsif ($line =~ /^# Improper Coeffs/) { 
        if (defined($next_line = <$FH_TOPO>)) {
            if ($next_line =~ /^#$/) {
                print $OUT "\n";
                handle_improper_coefficients();
                print $OUT "\n";
            }
        }
    } else {
        print $OUT $line;
    }
}

# write config file
print $CONFIG "# Create by $0.";
print $CONFIG "units           real\n";
print $CONFIG "neigh_modify    delay 2 every 1\n\n";
print $CONFIG "atom_style      full\n";
print $CONFIG "bond_style      harmonic\n";
print $CONFIG "angle_style     charmm\n";
print $CONFIG "dihedral_style  charmm\n";
print $CONFIG "improper_style  harmonic\n\n";
print $CONFIG "pair_style      lj/charmm/coul/long 8 10\n";
print $CONFIG "pair_modify     mix arithmetic\n";
print $CONFIG "kspace_style    pppm 1e-4\n\n";
print $CONFIG "read_data       ${output}.data\n\n";
print $CONFIG "special_bonds   charmm\n";
print $CONFIG "fix             1 all nve\n";
print $CONFIG "fix             2 all shake 1e-6 500 0 m 1.0\n";
print $CONFIG "velocity        all create 0.0 12345678 dist uniform\n\n";
print $CONFIG "thermo          1\n";
print $CONFIG "thermo_style    multi\n";
print $CONFIG "timestep        0.5\n\n";
print $CONFIG "dump            1 all atom 10 ${output}.dump\n";
print $CONFIG "dump_modify     1 image yes scale yes\n\n";
print $CONFIG "run             20\n";

# close files
close($NOTFOUND);
close($CONFIG);
close($OUT);
close($FH_FF);
close($FH_TOPO);

__END__

=head1 NAME

topo-ff - Parametrize a LAMMPS data file with a CHARMM force field

=head1 SYNOPSIS

=for pod2usage: 

topo-ff --topology TOP --ff FF --output OUT [--verbose] [--help] [--man]

=head1 OPTIONS

=over 

=item B<--topology>
    TOP - LAMMPS data file generated by Topotools

=item B<--ff>
    FF - CHARMM force field parameter file

=item B<--output>
    OUT - output file

=item B<--verbose>
    Print more debug information.

=item B<--help>
    Print a brief help message and exit.

=item B<--man>
    Print the manual page and exit.

=item B<--usage>
    Print the usage message and exit.


=back

=head1 DESCRIPTION

Parametrize the TOPOLOGY file for LAMMPS generated with 
the VMD plugin 'topotools' with a given CHARMM force field,
i.e. use the provided CHARMM topology (and parameter) file,
then write the completed LAMMPS data to an OUTPUT file.

=head1 TODOS

=over

=item a) dihedral weighting (last parameter, aka LJ 1-4 scaling
this is set to 0, but needs to be computed accordingly):
L<Parameters|http://lammps.sandia.gov/doc/dihedral_charmm.html>

=item b) take care of NBFIX terms

=item c) TIP3P for LAMMPS uses special LJ potentials:
L<Parameters|http://lammps.sandia.gov/doc/Section_howto.html#howto-7>

=item d) correct bounding box coordinates given by 'topotools'

=item e) update masses from CHARMM parameter file ('topotools' provides 
us already with the correct masses)) instead of 'topotools' provided
values

=item f) perl test cases

=back

=head1 NOTES

=over

=item Parameters in CHARMM file cannot be separated by empty lines,
they have the meaning to separate the different sections.

=item When generating a PSF/PDB file pair with VMD for instance,
the provided topology file (CHARMM) should match the 
provided parameter file (CHARMM) when using topo-ff.

=item When using 'topotools' then you need to retypebonds,
retypeangles, retypedihedrals and retypeimpropers,
then write the LAMMPS data with writelammpsdata.

=item LAMMPS currently cannot use CMAP corrections from FFs,
i.e. no cross terms are available and thus are ignored here.

=back

=head1 REQUIREMENTS              

Perl v5.000

=head1 BUGS

--

=head1 AUTHOR

Stephan Grein (grein@temple.edu)

=head1 ORGANIZATION

Temple University

=head1 REVISION

=over

=item Version 1.0

=item For more details cf. `git rev-parse HEAD`.

=back

=cut
