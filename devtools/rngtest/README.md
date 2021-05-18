Dieharder Results
=================

The results from the Dieharder test suite for the various RNGs in cpptraj are in the results.X files. The following tests were run:
```
 Test Number	                     Test Name	              Test Reliability
===============================================================================
  -d 0  	                  Diehard Birthdays Test	      Good
  -d 1  	                     Diehard OPERM5 Test	      Good
  -d 2  	          Diehard 32x32 Binary Rank Test	      Good
  -d 3  	            Diehard 6x8 Binary Rank Test	      Good
  -d 4  	                  Diehard Bitstream Test	      Good
  -d 8  	      Diehard Count the 1s (stream) Test	      Good
  -d 9  	        Diehard Count the 1s Test (byte)	      Good
  -d 10  	                Diehard Parking Lot Test	      Good
  -d 11  	Diehard Minimum Distance (2d Circle) Test	      Good
  -d 12  	Diehard 3d Sphere (Minimum Distance) Test	      Good
  -d 13  	                    Diehard Squeeze Test	      Good
  -d 15  	                       Diehard Runs Test	      Good
  -d 16  	                      Diehard Craps Test	      Good
  -d 17  	            Marsaglia and Tsang GCD Test	      Good
  -d 100  	                        STS Monobit Test	      Good
  -d 101  	                           STS Runs Test	      Good
  -d 102  	           STS Serial Test (Generalized)	      Good
  -d 201  	   RGB Generalized Minimum Distance Test	      Good
  -d 202  	                   RGB Permutations Test	      Good
  -d 203  	                     RGB Lagged Sum Test	      Good
  -d 204  	        RGB Kolmogorov-Smirnov Test Test	      Good
  -d 205  	                       Byte Distribution	      Good
  -d 206  	                                 DAB DCT	      Good
  -d 207  	                      DAB Fill Tree Test	      Good
  -d 208  	                    DAB Fill Tree 2 Test	      Good
  -d 209  	                      DAB Monobit 2 Test	      Good
```

- 0: Marsaglia: 44 failed
- 1: C stdlib: 37 failed
- 2: Mersenne Twister (mt19937 from C++ std): 2 failed.
- 3: Permuted congruential generator 32: 2 failed.
- 4: Xoshiro128++: 2 failed.

