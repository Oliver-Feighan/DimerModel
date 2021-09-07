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
Mg 15.77432 1.41642 1.19144
C 15.54581 3.85904 -1.42780
C 15.28739 -0.97056 -1.19600
C 15.52078 -0.72160 3.52210
C 15.10248 3.98170 3.45277
N 15.46931 1.55123 -0.99512
C 15.46319 2.52865 -1.88897
C 15.44370 2.04990 -3.36544
C 14.79209 0.59966 -3.13940
C 15.23605 0.35406 -1.69883
C 13.27292 0.61016 -3.17582
C 16.88421 1.94320 -3.97827
C 16.92558 1.78127 -5.49940
H 18.01351 2.41675 -6.23767
N 15.93966 -0.56070 1.14409
C 15.72765 -1.36987 0.08504
C 15.95957 -2.73863 0.57174
C 16.16479 -2.68426 2.04471
C 15.85563 -1.29216 2.30272
C 16.04958 -3.92390 -0.34502
C 16.55986 -3.71759 3.05582
O 16.44606 -3.50579 4.31480
C 16.99581 -5.01855 2.52768
N 15.16271 1.55712 3.17348
C 15.19773 0.48597 3.98547
C 14.97754 0.88534 5.42506
C 15.02825 2.48371 5.44007
C 15.18697 2.74608 3.92634
C 13.66713 0.41207 6.05374
C 16.03939 3.14663 6.41515
C 17.45397 2.85311 5.96996
N 15.57294 3.56523 1.09103
C 15.26452 4.44359 2.15285
C 15.12229 5.79479 1.57281
C 15.33683 5.60348 0.19954
C 15.57461 4.22683 -0.06082
C 14.78501 6.97688 2.29046
C 15.40571 6.25144 -1.09046
O 15.37246 7.41822 -1.47938
C 15.41890 5.13163 -2.21839
C 16.51701 5.46479 -3.19085
O 17.70769 5.43906 -2.83276
O 15.99881 5.54547 -4.42993
C 17.05138 5.86316 -5.43142
H 14.91616 -1.77940 -1.82597
H 15.52889 -1.40016 4.37825
H 14.99040 4.85045 4.10499
H 14.81909 2.63932 -4.02590
H 15.23142 -0.13885 -3.80495
H 12.87814 1.61382 -3.23017
H 12.95341 0.12691 -2.26180
H 12.91106 -0.02079 -3.99270
H 17.33233 0.98811 -3.67540
H 17.48500 2.78303 -3.68780
H 15.95948 2.01139 -5.96484
H 17.09051 0.72636 -5.73172
H 15.69335 -3.68033 -1.35357
H 15.33491 -4.66792 0.02316
H 17.04662 -4.36924 -0.36106
H 16.19939 -5.44898 1.92527
H 17.32158 -5.59264 3.39726
H 17.84601 -4.85625 1.85579
H 15.76885 0.50269 6.05791
H 14.06223 2.91927 5.73490
H 13.81517 -0.19212 6.96868
H 13.11609 -0.15031 5.31170
H 12.96588 1.20405 6.31109
H 15.93454 2.72135 7.41630
H 15.91131 4.22008 6.50622
H 17.64382 1.86065 6.38957
H 18.14377 3.59116 6.37433
H 17.51319 2.76650 4.87879
H 15.43821 7.75892 1.92683
H 14.89285 6.83479 3.36375
H 13.73886 7.09273 1.99480
H 14.43704 5.13341 -2.69329
H 16.72514 6.76207 -5.96360
H 17.18133 4.94236 -6.03057
H 18.01922 6.17757 -5.05875

