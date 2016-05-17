# topo-ff
Parametrize a LAMMPS data file with a CHARMM force field

## Usage

Create a PSF/PDB for a PDB file from the Protein Data Base.
Afterwards write out LAMMPS data file with the VMD plugin
'topotools'. Get the CHARMM22 parameter file, then invoke
topo-ff by:

```
./topo-ff.pl
```

This will print the help message:
```
topo-ff --topology TOP --ff FF --output OUT
```

Note, that you need to 
```chmod +x topo-ff.pl ```
previously or invoke by
``` perl topo-ff.pl ```.

You need to provide a topology file (TOP) and a
CHARMM22 force field (FF) and one output name for
the LAMMPS data file and input script (OUT).

There will be written two files with suffix *.data*
and *.in*.

## Requirements
Perl 5 is required only.

## Documentation
Invoke
```
./topo-ff.pl --man
```

## Help
Invoke
```
./topo-ff.pl --help
```

# (Short) usage
Invoke
```
./topo-ff.pl --usage
```

