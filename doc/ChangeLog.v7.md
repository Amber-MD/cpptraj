Changes in version 7.0 of CPPTRAJ. 2026
=======================================

The major change in this version is the introduction of the ability of CPPTRAJ to fully build and parameterize systems in the same way that LEaP does.

New Commands
============
- `build` - Allow full build and parameterization of systems for MD, including CMAP, LJ 12-6-4, solvation, and ions.

- `source` - Read (limited) leaprc files.

- `mutate` - Mutate residues from one kind to another, keeping only common atoms.

- `desc` - DEBUG: describe a selection in the same manner as leap.

New Keywords
============
- FIXME check this: Introduce the `byatom` keyword to `checkchirality` action to check the chirality of any atoms in specified mask that are chiral centers.

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

Add `gb` to `change` command to chang GB radii.

Add `complete` keyword to `zmatrix` command to get complete zmatrix instead of minimal zmatrix.

New functionality
=================
- Read topology/coordinates with the 'readdata' command.

- More complete read of Amber library/prep files.

- Read CMAP energy term from Amber MDOUT file.

- The `permutedihedrals` command will now check rings.

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
