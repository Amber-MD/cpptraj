Changes in version 7.0 of CPPTRAJ. 2026
=======================================

The major change in this version is the introduction of the ability of CPPTRAJ to fully build and parameterize systems in the same way that LEaP does.

New Commands
============
- `build` - Allow full build and parameterization of systems for MD, including CMAP, LJ 12-6-4, solvation, and ions. Includes `prepareforleap` and `fixatomorder` functionality. Can create macromolecular assemblies from PDBs with multiple models. Can solvate with a target # of solvent molecules.

- `source` - Read (limited) leaprc files.

- `mutate` - Mutate residues from one kind to another, keeping only common atoms.

- `desc` - DEBUG: describe a selection in the same manner as leap.

New functionality
=================
- Read topology/coordinates with the 'readdata' command.

- More complete read of Amber library/prep files.

- Read CMAP energy term from Amber MDOUT file.

- The `permutedihedrals` command will now check rings.

- The `graft` command will now use CONNECT atom information if it is present.

- The `energy` command will now calculate CMAP energies (if CMAP parameters are present) and/or LJ 12-6-4 contributions to VDW energy (if `lj1264` is specified).

- The `sequence` command now generates better geometries around bonds linking residues.

New Keywords
============
- Introduce ring intersection check to `checkstructure` action.
```
{noringcheck | [ringshortdist <rsdist>] [ringdcut <ringdcut>]
               [ringacut <ringacut>]}]
```
 Ring-bond intersections are detected when the ring center to bond center is less than <rsdist>, or less than <ringdcut> and the angle between the bond vector and ring perpendicular vector is less than <ringacut>.

- Ignore extra points by default in the `checkstructure action.
```
[{checkxp | xpmask <xpmask>}]
```
If `checkxp` is specified, extra points will be included. By default, extra points will be ignored as determined by 'xpmask'; the default mask to select extra points is '&@\\XP'

Add `namemap <mapset>` keywords for `change atomname from <mask>` to use an atom name map for determining new atom names.

Add `gb` to `change` command to change GB radii.

Add `complete` keyword to `zmatrix` command to get complete zmatrix instead of minimal zmatrix.

Add `atombondorder` keyword to mol2 output; sort bonds by first atom index (similar to how leap writes mol2 files).

Add `lj1264` keyword to the `energy` command to enable LJ 12-6-4 energy calculation.

New file formats recognized
===========================
- Read Amber parameter (force field) files.

- Read Amber force field modification (frcmod) files.

- Read leaprc files (limited).

- Read atom/pdb residue name maps.

Misc Fixes
==========
- Fix timing percentage of solute-solvent hydrogen bonds in `hbond` action.

- `replicatecell` action will correctly regenerate parameters as needed.

- Fix small memory leak when calculating modes.

- Ensure box information is properly updated if needed when using `crdaction`.

- Improve memory allocation in the Frame class to avoid frequent reallocation of memory.

- Improve handling of very large topology files.

- Better handling of very small boxes with pairlist code.

Not Yet Ready
=============
- Introduce the `byatom` keyword to `checkchirality` action to check the chirality of any atoms in specified mask that are chiral centers.

- Cannot yet build topology files for polarizable force fields (IPOL > 0).

ModXNA Compatibility
====================
The old `sequence` command functionality and some new features have been added to ensure backwards compatibility with ModXNA (V1.9.2 as of right now).

The following files/classes/functions have been retained from CPPTRAJ V6.31.0 with very few changes: `Structure/BuildAtom.h`, `Structure/Model.h`, `Structure/Model.cpp`

The `OldBuilder` class in `Structure/OldBuilder.cpp` and `Structure/OldBuilder.h` is the `Builder` class from CPPTRAJ V6.31.0.

The `Zmatrix::SetupICsAroundBond()` function has been retained from CPPTRAJ V6.31.0 `Zmatrix` class so that `sequence` can be made to work the way it used to.

Default bond Req parameters are set when reading in mol2 files. This is needed so that the RK/REQ columns are printed for the `bonds` command, which is expected/needed for modxna.sh to determine the name of the H atom bonded to O3' (`--3cap` flag).

ModXNA metadata (head/tail/anchor atom information) will be read from the title of Mol2 files if present and retained in the Topology. This metadata will be written back out to Mol2 files if another title has not been specified. This allows the `sequence` command to automatically recognized when units have ModXNA information and switch back to the old `sequence` mode automatically. 
