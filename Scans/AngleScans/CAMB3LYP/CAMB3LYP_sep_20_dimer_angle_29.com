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
Mg 10.51741 0.94752 0.80667
C 10.50224 -2.36152 2.19559
C 9.30777 2.04809 3.80001
C 10.09314 3.79611 -0.52486
C 10.62083 -0.52761 -2.34912
N 10.01189 -0.11549 2.67929
C 10.10552 -1.36426 3.11075
C 9.83156 -1.54271 4.62821
C 8.93078 -0.23458 4.86540
C 9.46905 0.64099 3.73565
C 7.44962 -0.47153 4.62176
C 11.14690 -1.49701 5.48249
C 10.99089 -1.97070 6.92922
H 12.09811 -2.69444 7.54787
N 10.27713 2.73216 1.64088
C 9.79445 2.99601 2.87323
C 9.79843 4.46150 3.00160
C 10.16782 5.04678 1.68414
C 10.17353 3.85560 0.85881
C 9.54981 5.16490 4.30423
C 10.45438 6.44853 1.23785
O 10.52162 6.75940 -0.00379
C 10.56131 7.46429 2.29534
N 10.16353 1.55004 -1.15141
C 10.07045 2.85363 -1.46749
C 10.09113 3.05490 -2.96397
C 10.46243 1.62836 -3.58423
C 10.50650 0.79193 -2.28661
C 8.78797 3.55634 -3.58582
C 11.68493 1.57352 -4.54087
C 12.95564 1.85280 -3.77101
N 10.74110 -1.05756 0.03359
C 10.73116 -1.45540 -1.32129
C 10.80094 -2.93091 -1.34231
C 10.82360 -3.29312 0.01297
C 10.75123 -2.12681 0.82167
C 10.78668 -3.74788 -2.50784
C 10.88156 -4.39749 0.94341
O 11.04110 -5.61305 0.83999
C 10.54851 -3.85062 2.39819
C 11.57984 -4.39688 3.34713
O 12.77271 -4.05993 3.24652
O 10.95792 -5.05145 4.34470
C 11.93868 -5.60208 5.31770
H 8.71633 2.46577 4.61527
H 10.05728 4.76043 -1.03698
H 10.75603 -1.05772 -3.29430
H 9.27030 -2.43256 4.88692
H 9.13891 0.21851 5.83119
H 7.26047 -1.45319 4.21347
H 7.14011 0.29483 3.92324
H 6.88283 -0.28963 5.53936
H 11.42451 -0.44947 5.65667
H 11.93151 -2.05006 5.00379
H 10.04610 -2.50448 7.08856
H 10.91508 -1.09288 7.57551
H 9.14310 4.48104 5.05943
H 8.74445 5.88648 4.12937
H 10.42974 5.69894 4.66920
H 9.63463 7.49295 2.86358
H 10.85717 8.38449 1.78781
H 11.34936 7.16081 2.99354
H 10.85318 3.77046 -3.24744
H 9.64078 1.22208 -4.19223
H 8.90960 4.49743 -4.15470
H 8.05876 3.68079 -2.79633
H 8.29154 2.85063 -4.24970
H 11.60544 2.35352 -5.30207
H 11.78512 0.62609 -5.05983
H 12.98683 2.94585 -3.73443
H 13.81880 1.45071 -4.29762
H 12.87841 1.49052 -2.73929
H 11.54030 -4.51045 -2.36278
H 10.97884 -3.16357 -3.40530
H 9.75934 -4.12014 -2.47254
H 9.54160 -4.18481 2.65123
H 11.74380 -6.67623 5.39459
H 11.81617 -5.00077 6.23811
H 12.98429 -5.59670 5.03287

