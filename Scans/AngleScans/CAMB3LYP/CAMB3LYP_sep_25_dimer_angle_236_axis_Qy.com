%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 0.00000 0.00000 0.00000
C -0.21700 -2.21800 2.80100
C -1.25500 2.44200 2.02600
C -0.24500 1.85000 -2.56400
C 0.14600 -2.82400 -2.03700
N -0.61700 -0.00900 2.12600
C -0.59400 -0.90800 3.10500
C -0.93600 -0.30900 4.50200
C -1.77300 0.96400 4.03200
C -1.15500 1.16800 2.66000
C -3.26300 0.68500 3.86300
C 0.34500 0.11000 5.29900
C 0.12400 0.39700 6.78700
H 1.16800 0.03300 7.73300
N -0.18500 1.96100 -0.15200
C -0.70300 2.80000 0.77700
C -0.62000 4.15600 0.19300
C -0.18700 4.01500 -1.22500
C -0.21400 2.58300 -1.38000
C -0.88700 5.40600 0.97400
C 0.18100 5.01300 -2.27700
O 0.30300 4.69000 -3.50800
C 0.31000 6.41400 -1.83700
N -0.26500 -0.41300 -2.02900
C -0.27900 0.56900 -2.93700
C -0.19600 0.02400 -4.33200
C 0.12400 -1.53300 -4.18000
C 0.08600 -1.63200 -2.63500
C -1.45100 0.20500 -5.19700
C 1.37300 -2.07000 -4.92300
C 2.62600 -1.49000 -4.30600
N 0.15000 -2.14400 0.30300
C 0.17300 -3.13600 -0.68800
C 0.17300 -4.43800 0.00400
C 0.13300 -4.11100 1.35800
C 0.08300 -2.69500 1.50200
C 0.14900 -5.71500 -0.62300
C 0.09600 -4.62900 2.71600
O 0.20800 -5.74700 3.22100
C -0.25500 -3.42800 3.70100
C 0.71700 -3.48200 4.84800
O 1.93100 -3.25500 4.66900
O 0.03800 -3.55600 6.00300
C 0.94300 -3.58300 7.18300
H -1.85900 3.22000 2.49700
H -0.21800 2.45300 -3.47400
H 0.27600 -3.74700 -2.60500
H -1.55400 -0.94900 5.13200
H -1.58500 1.83100 4.66500
H -3.48400 -0.37700 3.96900
H -3.51100 1.02700 2.85800
H -3.85200 1.30100 4.54300
H 0.66900 1.09600 4.96600
H 1.13000 -0.63900 5.18900
H -0.84400 0.02600 7.12500
H 0.06200 1.47700 6.92100
H -1.34400 5.18600 1.93800
H -1.65100 5.96300 0.43200
H 0.00700 6.02000 1.08400
H -0.63500 6.74400 -1.40500
H 0.66500 6.97000 -2.70400
H 1.04800 6.46600 -1.03600
H 0.60800 0.50100 -4.89100
H -0.69200 -2.15600 -4.54700
H -1.26700 0.74000 -6.12900
H -2.20400 0.71400 -4.59600
H -1.96700 -0.72000 -5.45300
H 1.35600 -1.75900 -5.96700
H 1.45200 -3.15700 -4.90200
H 2.71400 -0.52200 -4.79800
H 3.48700 -2.12100 -4.52900
H 2.49500 -1.30300 -3.24000
H 0.86700 -6.33300 -0.08400
H 0.40400 -5.64200 -1.68000
H -0.89100 -5.99900 -0.46300
H -1.28300 -3.58100 4.02800
H 0.69800 -4.47800 7.75500
H 0.83300 -2.61900 7.68000
H 2.00000 -3.75400 6.97600
Mg 13.12844 1.18575 0.99915
C 11.45035 -1.62932 -0.46321
C 11.64375 3.12008 -1.38595
C 14.93483 3.60392 1.97732
C 15.31826 -1.07410 2.49697
N 11.74556 0.69339 -0.65539
C 11.12834 -0.40042 -1.07605
C 10.04922 -0.15199 -2.16386
C 10.61353 1.23224 -2.75067
C 11.35368 1.73823 -1.51421
C 11.62500 1.04397 -3.86907
C 8.61472 0.02526 -1.55328
C 7.47144 -0.03727 -2.56854
H 6.20964 -0.64778 -2.15913
N 12.91862 3.13967 0.72208
C 12.26517 3.75290 -0.28699
C 12.41873 5.19622 -0.04737
C 13.35504 5.38524 1.09381
C 13.79837 4.02102 1.29980
C 11.65094 6.23439 -0.81296
C 13.79038 6.59560 1.86304
O 14.79343 6.56192 2.66028
C 13.06621 7.84452 1.58506
N 14.99688 1.27293 1.90670
C 15.53168 2.44993 2.27629
C 16.74510 2.24677 3.15192
C 16.77180 0.68463 3.49341
C 15.56714 0.21974 2.64603
C 18.08401 2.65637 2.53879
C 16.81351 0.29332 4.99601
C 15.51315 0.67609 5.66530
N 13.26048 -0.96290 1.18328
C 14.27362 -1.70280 1.83150
C 13.99013 -3.13257 1.59103
C 12.83552 -3.12457 0.79398
C 12.43972 -1.78397 0.53778
C 14.78144 -4.22531 2.04476
C 11.86762 -3.94517 0.10224
O 11.64333 -5.15296 0.03394
C 11.00029 -3.01144 -0.84764
C 9.55113 -3.35866 -0.64297
O 8.99829 -3.14488 0.45032
O 9.00352 -3.68110 -1.82894
C 7.56380 -4.02262 -1.67934
H 11.41215 3.77874 -2.22345
H 15.53260 4.39877 2.42939
H 15.89811 -1.84319 3.01169
H 10.01254 -0.90020 -2.94656
H 9.80452 1.90858 -3.01415
H 11.88664 0.00512 -4.00500
H 12.49494 1.61836 -3.57870
H 11.25622 1.50114 -4.79175
H 8.49951 1.05943 -1.20404
H 8.44277 -0.68828 -0.77103
H 7.79899 -0.44229 -3.53374
H 7.15601 0.98317 -2.79923
H 11.16604 5.80463 -1.69816
H 12.38431 6.93886 -1.22010
H 10.94415 6.78271 -0.18650
H 13.15709 8.08727 0.52900
H 13.46009 8.57483 2.29455
H 12.00032 7.68169 1.77996
H 16.65429 2.80668 4.07461
H 17.66053 0.19092 3.07359
H 18.63145 3.40280 3.14477
H 17.90123 3.03584 1.54223
H 18.77576 1.83702 2.35074
H 17.60223 0.84812 5.51008
H 17.00130 -0.76215 5.16243
H 15.65836 1.73500 5.89922
H 15.36030 0.08959 6.56901
H 14.67378 0.60557 4.96377
H 14.09457 -4.97344 2.41760
H 15.48644 -3.91382 2.81268
H 15.28395 -4.50192 1.11404
H 11.32451 -3.19535 -1.87265
H 7.42667 -5.02187 -2.10406
H 7.00882 -3.19247 -2.15535
H 7.18453 -4.16611 -0.67435

