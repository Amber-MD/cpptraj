TestU01 Results
===============

These are results from the TestU01 test suite for the various RNGs in cpptraj.

Marsaglia
=========
The Marsaglia generator didnt complete tests.

C Stdlib
========
```
 Version:          TestU01 1.2.3
 Generator:        CPPTRAJ
 Number of statistics:  144
 Total CPU time:   00:32:50.33
 The following tests gave p-values outside [0.001, 0.9990]:
 (eps  means a value < 1.0e-300):
 (eps1 means a value < 1.0e-15):

       Test                          p-value
 ----------------------------------------------
  1  SerialOver, t = 2                eps  
  2  SerialOver, t = 4                eps  
  3  CollisionOver, t = 2             eps  
  5  CollisionOver, t = 4             eps  
  7  CollisionOver, t = 8             eps  
  9  CollisionOver, t = 20            eps  
 11  BirthdaySpacings, t = 2          eps  
 12  BirthdaySpacings, t = 3          eps  
 13  BirthdaySpacings, t = 4          eps  
 14  BirthdaySpacings, t = 7          eps  
 18  ClosePairs NP, t = 2            2.0e-4
 18  ClosePairs mNP, t = 2          4.6e-66
 18  ClosePairs mNP1, t = 2         8.4e-70
 18  ClosePairs NJumps, t = 2         eps  
 19  ClosePairs NP, t = 3            1.2e-7
 19  ClosePairs mNP, t = 3         1.2e-134
 19  ClosePairs mNP1, t = 3        1.3e-142
 19  ClosePairs NJumps, t = 3         eps  
 20  ClosePairs NP, t = 7           7.7e-10
 20  ClosePairs mNP, t = 7          1.8e-79
 20  ClosePairs mNP1, t = 7        1.9e-240
 20  ClosePairs NJumps, t = 7         eps  
 23  SimpPoker, d = 16                eps  
 25  SimpPoker, d = 64                eps  
 26  SimpPoker, d = 64                eps  
 27  CouponCollector, d = 4           eps  
 29  CouponCollector, d = 16          eps  
 30  CouponCollector, d = 16         6.2e-8
 31  Gap, r = 0                       eps  
 32  Gap, r = 27                      eps  
 33  Gap, r = 0                       eps  
 34  Gap, r = 22                      eps  
 41  MaxOft, t = 5                    eps  
 41  MaxOft AD, t = 5              3.2e-157
 42  MaxOft, t = 10                   eps  
 42  MaxOft AD, t = 10              1.8e-79
 43  MaxOft, t = 20                   eps  
 43  MaxOft AD, t = 20              1 - eps1
 44  MaxOft, t = 30                   eps  
 44  MaxOft AD, t = 30              1 - eps1
 45  SampleProd, t = 10             1 - eps1
 46  SampleProd, t = 30             1 - eps1
 47  SampleMean                       eps  
 48  SampleCorr                     1 - eps1
 49  AppearanceSpacings, r = 0      1 - eps1
 51  WeightDistrib, r = 0             eps  
 52  WeightDistrib, r = 8             eps  
 53  WeightDistrib, r = 16            eps  
 54  WeightDistrib, r = 24            eps  
 55  SumCollector                     eps  
 56  MatrixRank, 60 x 60              eps  
 58  MatrixRank, 300 x 300            eps  
 60  MatrixRank, 1200 x 1200          eps  
 62  Savir2                           eps  
 65  RandomWalk1 H (L = 90)           eps  
 65  RandomWalk1 M (L = 90)           eps  
 65  RandomWalk1 J (L = 90)           eps  
 65  RandomWalk1 R (L = 90)           eps  
 65  RandomWalk1 C (L = 90)           eps  
 67  RandomWalk1 H (L = 1000)         eps  
 67  RandomWalk1 M (L = 1000)         eps  
 67  RandomWalk1 J (L = 1000)         eps  
 67  RandomWalk1 R (L = 1000)         eps  
 67  RandomWalk1 C (L = 1000)         eps  
 69  RandomWalk1 H (L = 10000)        eps  
 69  RandomWalk1 M (L = 10000)        eps  
 69  RandomWalk1 J (L = 10000)        eps  
 69  RandomWalk1 R (L = 10000)        eps  
 69  RandomWalk1 C (L = 10000)        eps  
 71  LinearComp, r = 0              1 - eps1
 71  LinearComp, r = 0              1 - eps1
 73  LempelZiv                      1 - eps1
 74  Fourier3, r = 0                  eps  
 76  LongestHeadRun, r = 0            eps  
 76  LongestHeadRun, r = 0          1 -  8.1e-9
 80  HammingWeight2, r = 0            eps  
 82  HammingCorr, L = 30              eps  
 83  HammingCorr, L = 300             eps  
 84  HammingCorr, L = 1200            eps  
 85  HammingIndep, L = 30             eps  
 87  HammingIndep, L = 300            eps  
 89  HammingIndep, L = 1200           eps  
 91  Run of bits, r = 0               eps  
 95  AutoCor, d = 30                1 - eps1
 ----------------------------------------------
 All other tests were passed
```

Mersenne Twister
================
```
 Version:          TestU01 1.2.3
 Generator:        CPPTRAJ
 Number of statistics:  144
 Total CPU time:   00:31:10.94
 The following tests gave p-values outside [0.001, 0.9990]:
 (eps  means a value < 1.0e-300):
 (eps1 means a value < 1.0e-15):

       Test                          p-value
 ----------------------------------------------
 67  RandomWalk1 C (L = 1000)        6.7e-4
 71  LinearComp, r = 0              1 - eps1
 72  LinearComp, r = 29             1 - eps1
 ----------------------------------------------
 All other tests were passed
```

PCG32
=====
```
 Version:          TestU01 1.2.3
 Generator:        CPPTRAJ
 Number of statistics:  144
 Total CPU time:   00:27:29.09
 The following tests gave p-values outside [0.001, 0.9990]:
 (eps  means a value < 1.0e-300):
 (eps1 means a value < 1.0e-15):

       Test                          p-value
 ----------------------------------------------
 44  MaxOft, t = 30                  3.6e-4
 ----------------------------------------------
 All other tests were passed
```

Xoshiro128++
============
```
 Version:          TestU01 1.2.3
 Generator:        CPPTRAJ
 Number of statistics:  144
 Total CPU time:   00:28:41.95

 All tests were passed
```
