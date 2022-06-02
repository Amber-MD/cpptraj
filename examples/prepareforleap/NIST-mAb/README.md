CPPTRAJ 'prepareforleap' command example with NIST-mAb
======================================================

An example of how to use the 'prepareforleap' command to quickly parameterize glycoforms of the NIST monoclonal antibody.

https://nvlpubs.nist.gov/nistpubs/jres/126/jres.126.012.pdf

mAb8671_2021_0526.pdb : The NIST-mAb structure.

RunExample.sh : The example script.

This will generate:

leap.ff.in : File containing LEaP input for creating bonds.

Prepared.pdb : PDB ready for processing by LEaP.

Prepared.mol2 : Mol2 of prepared PDB.

If LEaP is available, will also generate:

Mol.parm7 : Topology from LEaP.

Mol.rst7 : Coordinates from LEaP.
