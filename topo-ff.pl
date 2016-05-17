#!/usr/bin/env perl
## topo-ff - Parametrize a LAMMPS data file with a CHARMM force field
# TODO test cases and code cleanup

use 5.000;
use utf8;
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Copy qw(copy);
use constant VERSION => "0.1";

# CLI
my $topo;
my $ff; 
my $output;
my $verbose = 0;
my $no_dihedrals = 0;
my $no_angles = 0;
my $no_impropers = 0;
my $no_pairs = 0;
my $no_bonds = 0;
my $not_found = "not_found.txt";

GetOptions(
    'topology=s'        => \$topo,
    'ff=s'              => \$ff,
    'output=s'          => \$output,
    'verbose'           => \$verbose,
    'no-angles'         => \$no_angles,
    'no-bonds'          => \$no_bonds,
    'no-dihedrals'      => \$no_dihedrals,
    'no-pairs'          => \$no_pairs,
    'no-impropers'      => \$no_impropers,
    'usage'             => sub { pod2usage(2) },
    'help'              => sub { pod2usage(1) },
    'man'               => sub { pod2usage(-existatus => 0, -verbose => 2) },
);

pod2usage(2) unless defined($topo) and defined($ff) and defined($output);

# open files
open(my $FH_TOPO, "<:encoding(UTF-8)", "$topo")
  or die "Could not open toplogy file '$topo': $!";
open(my $FH_FF, "<:encoding(UTF-8)", "$ff")
  or die "Could not open force field file '$ff': $!";
open (my $OUT, ">:encoding(UTF-8)", "${output}.tmp")
    or die "Could not open output file '$output': $!";
open (my $CONFIG, ">:encoding(UTF-8)", "${output}.in")
    or die "Could not open config file '${output}.in': $!"; 
open (my $NOTFOUND, ">:encoding(UTF-8)", "$not_found")
    or die "Could not open file '$not_found': $!";

# coefficients from CHARMM parameter file
my @pairs;
my @bonds;
my @angles;
my @dihedrals;
my @impropers;

# handle pair coefficients
sub handle_pair_coefficients {
    print $NOTFOUND " Pair Coeffs\n\n";
    print $OUT " Pair Coeffs\n\n";
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
                    if ($#output+1 == 2) {
                       # print $OUT "$index $val1 $val2 $val1 $val2 # $line\n";
                        print $OUT "$index $val1 $val2 # $line\n";
                    } else {
                        my $val3 = abs($output[2]);
                        my $val4 = 2.0*($output[2]*2.0**(-1/6));
                        print $OUT "$index $val1 $val2 $val3 $val4 # $line\n";
                    }
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
    print $NOTFOUND " Bond Coeffs\n\n";
    print $OUT " Bond Coeffs\n\n";
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
    print $OUT " Angle Coeffs\n\n";
    print $NOTFOUND " Angle Coeffs\n\n";
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
    print $NOTFOUND " Dihedral Coeffs\n\n";
    print $OUT " Dihedral Coeffs\n\n";
    my $weighting = 1; # default LJ 1-4 weighting 1, in rings different value!
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
                    print $OUT "$index $output[0] " . int($output[1]) . " " .  int($output[2]) . " $weighting # $line \n";
                    $index++;
                    $found = 1;
                    last; 
                }

                 if ($dihedral =~ /$fromto[0]\s*$fromto[1]\s*$fromto[2]\s*$fromto[3]/) {
                    my @output = ($dihedral=~ /(\d+\.\d+)\s*(\d+)\s*(\d+\.\d+)/);
                    print $OUT "$index $output[0] " . int($output[1]) . " " .  int($output[2]) . " $weighting # $line \n";
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
                        print $OUT "$index $output[0] " . int($output[1]) . " " .  int($output[2]) . " $weighting # $line \n";
                        $index++;
                        $found = 1;
                        last; 
                    }

                    if ($dihedral =~ /X\s*$fromto[2]\s*$fromto[1]\s*X/) {
                        my @output = ($dihedral=~ /(\d+\.\d+)\s*(\d+)\s*(\d+\.\d+)/);
                        print $OUT "$index $output[0] " . int($output[1]) . " " .  int($output[2]) . " $weighting # $line \n";
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
    print $OUT " Improper Coeffs\n\n";
    print $NOTFOUND " Improper Coeffs\n\n";
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
                     print $OUT "$index $output[0] $output[2] # $line \n";
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
        } else {
            return;
        }
    }
}

# read coefficients for different sections in CHARMM parameter file
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

# read in all coefficients for each section specified in CHARMM parameter file
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
    print "Information from CHARMM parameter file:\n";
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
    print "****************************\n\n\n";
}

# write data file header
my $header = <$FH_TOPO>;
chomp($header);
$header .= " (parametrized with CHARMM force field parameters by $0 v" . VERSION . ")\n\n";
print $OUT $header;

# write data file
while (my $line = <$FH_TOPO>) {
    my $next_line;
    if ($line =~ /^# Pair Coeffs/ && ! $no_pairs) {
        if (defined ($next_line = <$FH_TOPO>)) {
           if ($next_line =~ /^#$/) {
                print $OUT "\n";
                handle_pair_coefficients();
                print $OUT "\n";
           }
        }
    } elsif ($line =~ /^# Bond Coeffs/ && ! $no_bonds) {
        if (defined($next_line = <$FH_TOPO>)) {
            if ($next_line =~ /^#$/) {
                print $OUT "\n";
                handle_bond_coefficients();
                print $OUT "\n";
            }
        }
    } elsif ($line =~ /^# Angle Coeffs/ && ! $no_angles) { 
        if (defined($next_line = <$FH_TOPO>)) {
            if ($next_line =~ /^#$/) {
                print $OUT "\n";
                handle_angle_coefficients();
                print $OUT "\n";
            }
        }
    } elsif ($line =~ /^# Dihedral Coeffs/ && ! $no_dihedrals) { 
        if (defined($next_line = <$FH_TOPO>)) {
            if ($next_line =~ /^#$/) {
                print $OUT "\n";
                handle_dihedral_coefficients();
                print $OUT "\n";
            }
        }
    } elsif ($line =~ /^# Improper Coeffs/ && ! $no_impropers) { 
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
print $CONFIG "# Generated by $0.\n\n";
print $CONFIG "units           real\n";
print $CONFIG "neigh_modify    delay 2 every 1\n\n";
print $CONFIG "atom_style      full\n";
print $CONFIG "bond_style      " . (!$no_bonds ? "harmonic" : "none") . "\n";
print $CONFIG "angle_style     " . (!$no_angles ? "charmm" : "none") . "\n";
print $CONFIG "dihedral_style  " . (!$no_dihedrals ? "charmm" : "none") . "\n";
print $CONFIG "improper_style  " . (!$no_impropers ? "harmonic" : "none") . "\n";
print $CONFIG "pair_style      " . (!$no_pairs ? "lj/charmm/coul/long 8 10" : "none") . "\n";
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

sub read_entities {
    my $entity = shift;
    my $name = shift;
    my $IN = shift;
    while (my $line = <$IN>) {
        if ($line =~ /^\s*$name/) {
            my $next_line = <$IN>;
            if ($next_line =~ /^\s*$/) { # empty line marks start
                while (my $very_next_line = <$IN>) {
                    if ($very_next_line =~ /^\s*$/) { # empty line marks end
                        last;
                    } else {
                        push @$entity, $very_next_line;
                    }
                }
            }
            last;
        }
    }
    seek $IN, 0, 0;
}

sub print_entities_info {
    my $entity = shift;
    my $name = shift;
    print "\t$name: ";
    print scalar @$entity . "\t \n";
}

# correct dihedral weighting parameter
sub correct_dihedral_screening {
    # open input and output files
    open (my $IN, "<:encoding(UTF-8)", "${output}.tmp")
    or die "Could not open output file '$output': $!";

    open (my $OUT, ">:encoding(UTF-8)", "${output}.data")
    or die "Could not open output file '$output': $!";

    # sections in LAMMPS data file
    my @dihedrals;
    my @dihedral_coeffs;
    my @angles;
    my @angle_coeffs;
    my @atoms;

    my %types = ( "Dihedrals"       => \@dihedrals,
                  "Dihedral Coeffs" => \@dihedral_coeffs,
                  "Angles"          => \@angles,
                  "Angle Coeffs"    => \@angle_coeffs,
                  "Atoms"           => \@atoms,
                );

    # internal statistics (angle and dihedrals w coefficients during correction)
    if ($verbose) {
        print "Information from LAMMPS data file\n";
        print "(during correcting dihedral coefficients):\n";
        print "*******************************************\n";
    }

    for my $type (keys %types) {
        read_entities($types{$type}, $type, $IN);
        print_entities_info($types{$type}, $type) unless !$verbose;
    }

    if ($verbose) {
        print "*******************************************\n\n";
    }

    # correct the dihedral coefficients
    my %hash;
    my $hash_id;
    my $id1;
    my $id2;
    my $first;
    my $last;
    my @ids;

    for my $dihedral (@dihedrals) {
        my @columns = split " ", $dihedral;
        $id1 = $columns[1]; # dihedral_coefficient type
        $first = $columns[2]; # first atom defining dihedral
        $last  = $columns[5]; # last atom defining dihedral
        ($first, $last) = ($last, $first) if ($first>$last); # swap

        if (!defined($id2 = $hash{$hash_id = $first . " " . $last})) {
           $hash{$hash_id} = $id1;
         } else {
            push @ids, $id1-1;
            push @ids, $id2-1;
        }
    }

    ## shared 1-4 in 6-membered rings
    # n 6-membered rings, the same 1-4 interaction would be computed
    # twice (once for the clockwise 1-4 pair in dihedral 1-2-3-4 and
    # once in the counterclockwise dihedral 1-6-5-4) and thus the
    # weighting factor has to be 0.5 in this case.
    my @unique_ids = do { my %seen; grep { !$seen{$_}++ } @ids };
    for my $id (@unique_ids) {
        $dihedral_coeffs[$id] =~ s/(.*)\d(\s*#\s*.*)/${1}0.5${2}/;
    }

    @ids = ();
    for my $angle (@angles) {
        my @columns = split " ", $angle;
        $id1 = $columns[1]; # dihedral_coefficient type
        $first = $columns[2]; # first atom defining dihedral
        $last = $columns[4]; # first atom defining dihedral
        ($first, $last) = ($last, $first) if ($first>$last); # swap
        if (defined(($id1 = $hash{$first." ".$last}))) {
            push @ids, $id1-1;
        }
    }

    ## non-shared 1-4 in 5-membered rings
    # In 4-membered or 5-membered rings, the 1-4 dihedral also is
    # counted as a 1-2 or 1-3 interaction when going around the
    # ring in the opposite direction and thus the weighting factor
    # is 0.0, as the 1-2 and 1-3 exclusions take precedence.
    @unique_ids = do { my %seen; grep { !$seen{$_}++ } @ids };
    for my $id (@unique_ids) {
        $dihedral_coeffs[$id] =~ s/(.*)\d(\s*#\s*.*)/${1}0${2}/;
    }
    
    # write final corrected LAMMPS data file
    seek $IN, 0, 0;
    while (my $line = <$IN>) {
        if ($line =~ /^\s*Dihedral Coeffs/) {
            print $OUT "Dihedral Coeffs\n\n";
            my $next_line = <$IN>;
            my $counter = 0;
            if ($next_line =~ /^\s*$/) { # empty line marks start
                while (my $very_next_line = <$IN>) {
                     if ($very_next_line =~ /^\s*$/) { # empty line marks end
                         last;
                      } else {
                         print $OUT $dihedral_coeffs[$counter];
                         $counter++;
                     }
                 }
            }
        } else {
           print $OUT $line;
        }
    }
    close($OUT);
    close($IN);
    
    return @atoms;
}

# correct bounding box coordinates
sub correct_bounding_box {
    # open input and output files
    open (my $IN, "<:encoding(UTF-8)", "${output}.data")
    or die "Could not open output file '$output': $!";

    open (my $OUT, ">:encoding(UTF-8)", "${output}.tmp")
    or die "Could not open output file '$output': $!";

    my $atoms = shift;
    my @x;
    my @y;
    my @z;
    for my $atom (@$atoms) {
        my @cols = split " ", $atom;
        push @x, $cols[4];
        push @y, $cols[5];
        push @z, $cols[6];
    }

    my @xsorted = sort { $a <=> $b } @x;
    my @ysorted = sort { $a <=> $b } @y;
    my @zsorted = sort { $a <=> $b } @z;

    if ($verbose) {
       print "Bounding Box dimensions:\n";
       print "*******************************************\n";
       print "x min: " . $xsorted[0] . " | x max: " . $xsorted[-1] . "\n";
       print "y min: " . $ysorted[0] . " | y may: " . $ysorted[-1] . "\n";
       print "z min: " . $zsorted[0] . " | z maz: " . $zsorted[-1] . "\n";
       print "*******************************************\n";
    }

    # write bounding box
    while (my $line = <$IN>) {
        if ($line =~ /^.*?xlo\s*xhi$/) {
            print $OUT " $xsorted[0] $xsorted[-1] xlo xhi\n";
        } elsif ($line =~ /^.*?ylo\s*yhi$/) {
            print $OUT " $ysorted[0] $ysorted[-1] ylo yhi\n";
        } elsif ($line =~ /^.*?zlo\s*zhi$/) {
            print $OUT " $zsorted[0] $zsorted[-1] zlo zhi\n";
        } else {
            print $OUT $line;
        }
    }

    close ($OUT);
    close ($IN);
    copy "${output}.tmp", "${output}.data";
}

my @atoms = correct_dihedral_screening;
correct_bounding_box \@atoms;


__END__

=head1 NAME

topo-ff - Parametrize a LAMMPS data file with a CHARMM force field

=head1 SYNOPSIS

=for pod2usage: 

topo-ff --topology TOP --ff FF --output OUT [--verbose] [--no-angles] [--no-dihedrals]
[--no-impropers] [--no-bonds] [--no-pairs] [--help] [--man] [--usage]

=head1 OPTIONS

=over 

=item B<--topology>
    TOP - LAMMPS data file generated by Topotools

=item B<--ff>
    FF - CHARMM force field parameter file

=item B<--output>
    OUT - output file basename

=item B<--verbose>
    Print more debug information.

=item B<--no-bonds>
    Do not write bond coefficients

=item B<--no-angles>
    Do not write angle coefficients

=item B<--no-dihedrals>
    Do not write dihedral coefficients

=item B<--no-impropers>
    Do not write improper coefficients

=item B<--no-pairs>
    Do not write pair coefficients

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
i.e. use the provided CHARMM parameter file, then write the 
completed LAMMPS data and run script to two OUTPUT files 
(with suffixes .data and .in).

=head1 NOTES

=over

=item A sane CHARMM parameter file is required, thus entries
within sections cannot be separated by newlines - newlines
have the meaning of a section separator within CHARMM.

=item When using the VMD plugin 'topotools' then one is required
to use retypebonds, retypeangles, retypedihedrals and retypeimpropers
commands to prepare the intermediate LAMMPS output with the appropriate
hints for the force field parametrization before using writelammpsdata.

=item Current development version of LAMMPS does not incldue cross-terms,
also known as CMAP corrections from the CHARMM force field, therefore,
topo-ff should be used with a CHARMM22 parameter file.

=item LAMMPS used LJ 1-4 interaction differently (specified
by the pair_coeff command but applied in the dihedral section)
The last parameter for the dihedral_coeff (weighting) is therefore
adjusted, cf. comments in the code and read the LAMMPS documentation
at L<Parameters|http://lammps.sandia.gov/doc/dihedral_charmm.html>

=item In CHARMM impropers are defined from the topology file. The VMD
plugin 'topotools' just guesses them for structures that look like they 
might be requiring them to stay planar. Make sure that the impropers
are defined the way you intended before using topo-ff.

=item There are overrides to the mixing rules for i, j non-bonded parameters
from the "NBFIX" section in the parameter stream files (e.g. CHARMM for ions).
Make sure that the parametrized LAMMPS data file is defined in the way you
intended before running serious molecular dynamics simulations with LAMMPS.

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

=item Version 0.1

=item For more details cf. `git rev-parse HEAD`.

=back

=cut
