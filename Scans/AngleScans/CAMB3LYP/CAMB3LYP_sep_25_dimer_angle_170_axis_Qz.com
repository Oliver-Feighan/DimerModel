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
Mg 13.14666 1.17826 0.99297
C 13.02401 2.89304 -2.15720
C 12.65581 -1.72098 -0.73563
C 12.79991 -0.30587 3.77587
C 12.48567 4.23732 2.52625
N 12.88858 0.76785 -1.16509
C 12.92168 1.49115 -2.27418
C 12.92138 0.65934 -3.58477
C 12.23392 -0.67786 -3.02072
C 12.64344 -0.56317 -1.55382
C 10.71637 -0.65233 -3.09636
C 14.37122 0.37986 -4.11568
C 14.43954 -0.15732 -5.54694
H 15.55565 0.25613 -6.39301
N 13.26979 -0.75048 1.44395
C 13.06148 -1.79488 0.61504
C 13.25370 -3.00236 1.43309
C 13.43048 -2.58531 2.85067
C 13.14663 -1.16808 2.74585
C 13.33623 -4.38025 0.84316
C 13.78258 -3.33973 4.09689
O 13.64821 -2.81855 5.26014
C 14.20058 -4.73819 3.92063
N 12.49864 1.81913 2.86157
C 12.49401 0.98418 3.91550
C 12.25376 1.73378 5.20421
C 12.33898 3.28426 4.82168
C 12.53370 3.15785 3.29475
C 10.92100 1.45360 5.89823
C 13.34457 4.15325 5.62562
C 14.76106 3.73514 5.30292
N 12.99421 3.23693 0.35526
C 12.68378 4.35740 1.15668
C 12.58269 5.52312 0.25486
C 12.82051 4.99162 -1.02164
C 13.03340 3.58992 -0.92472
C 12.25693 6.85221 0.64667
C 12.92934 5.29584 -2.43029
O 12.92934 6.32897 -3.09842
C 12.94074 3.92984 -3.24289
C 14.06514 3.99196 -4.24008
O 15.24755 4.03727 -3.85740
O 13.57377 3.76910 -5.47266
C 14.65288 3.80973 -6.49524
H 12.27974 -2.65539 -1.15318
H 12.77606 -0.74925 4.77403
H 12.37950 5.24307 2.93840
H 12.32313 1.07521 -4.38660
H 12.67031 -1.56611 -3.47010
H 10.34473 0.31221 -3.40891
H 10.36814 -0.88689 -2.09898
H 10.35730 -1.46131 -3.73897
H 14.79226 -0.47647 -3.57330
H 14.98422 1.25582 -4.02882
H 13.48821 -0.03515 -6.07891
H 14.58607 -1.23934 -5.50487
H 13.00569 -4.39047 -0.20280
H 12.59828 -4.99719 1.36726
H 14.32345 -4.83149 0.96340
H 13.40721 -5.29250 3.42489
H 14.49626 -5.08220 4.91366
H 15.06742 -4.76247 3.25084
H 13.02369 1.50854 5.93193
H 11.37696 3.79513 4.97456
H 11.03746 1.09464 6.93826
H 10.37283 0.73273 5.30629
H 10.23216 2.29593 5.93258
H 13.21042 3.99309 6.69823
H 13.23811 5.21736 5.44309
H 14.92079 2.87593 5.96120
H 15.45853 4.53956 5.52762
H 14.84025 3.37793 4.26958
H 12.93418 7.50810 0.11597
H 12.34010 6.98085 1.72385
H 11.21969 6.90742 0.30556
H 11.96888 3.82882 -3.72749
H 14.35705 4.55246 -7.24259
H 14.77473 2.76654 -6.84260
H 15.61968 4.19159 -6.18875

