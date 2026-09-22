# ff19SB amino-acid parent library and ModXNA parent fragments for timap.
#
# Baseline: Amber ff19SB amino19.lib (the 20 amino acids plus common
# protonation variants). cpptraj ships this file so a noncanonical residue
# can be walked into the same N, H, CA, HA, side-chain, C, O order Amber
# uses, with no parent required:
#
#   readdata ncaa.lib name ncaa
#   timap ncaa[ncaa] aaorder out ncaa.lib
#
# To map an analog onto a specific ff19SB residue (TI against cysteine, ...):
#
#   readdata $CPPTRAJHOME/dat/timap/amino19.lib name FF
#   readdata ncaa.lib name ncaa
#   timap ncaa[ncaa] template FF[CYS] out ncaa.lib
#
# ModXNA consensus parent mol2 fragments also live here (see PARENTS.txt):
#
#   parm $CPPTRAJHOME/dat/timap/DAA.mol2 name DAA
#   timap analog template DAA out analog.lib
#
# Source: AmberTools ff19SB amino19.lib; ModXNA parents from modna-ff.
# See the timap section of the manual.
