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
Mg 8.93520 0.80690 0.66844
C 11.33182 -1.63386 -0.41691
C 11.16819 3.20358 -0.28251
C 7.04459 2.96509 2.02314
C 7.42485 -1.71316 2.54297
N 10.95552 0.67691 -0.22244
C 11.73853 -0.30314 -0.64759
C 12.99212 0.16933 -1.43172
C 13.11322 1.65001 -0.82198
C 11.65052 1.88065 -0.44782
C 13.94641 1.71186 0.44739
C 12.74962 0.20439 -2.98155
C 14.01577 0.35728 -3.82728
H 14.08178 -0.33271 -5.11259
N 8.84584 2.77364 0.41841
C 9.83961 3.58720 0.00391
C 9.27469 4.94552 0.00569
C 7.91695 4.88808 0.61227
C 7.89087 3.51350 1.07045
C 9.98017 6.12100 -0.60586
C 6.82784 5.90594 0.76718
O 5.82356 5.70551 1.53795
C 7.01826 7.18222 0.06266
N 7.57834 0.70450 2.24031
C 6.84756 1.77642 2.59370
C 5.75013 1.39452 3.55816
C 5.74479 -0.20437 3.59166
C 6.95177 -0.48085 2.66854
C 5.89808 1.93261 4.98117
C 4.39727 -0.91037 3.27791
C 4.02215 -0.68670 1.83062
N 9.17271 -1.33030 0.87833
C 8.45282 -2.18812 1.73855
C 9.04926 -3.53361 1.61065
C 10.09782 -3.35896 0.69494
C 10.16390 -1.99641 0.29670
C 8.64331 -4.69500 2.32653
C 11.17053 -4.01530 -0.01725
O 11.54323 -5.18091 -0.14449
C 12.09725 -2.90250 -0.67243
C 12.35025 -3.29207 -2.10307
O 11.41706 -3.32157 -2.92459
O 13.67755 -3.35353 -2.31489
C 13.97891 -3.72561 -3.72292
H 11.88174 4.02635 -0.33527
H 6.29893 3.64635 2.43950
H 6.94886 -2.57261 3.01988
H 13.89522 -0.39260 -1.22576
H 13.43054 2.36370 -1.57784
H 14.22243 0.72794 0.79647
H 13.32875 2.21392 1.18052
H 14.81750 2.35565 0.29554
H 12.23547 1.13837 -3.24248
H 12.20379 -0.66183 -3.30148
H 14.92590 0.17143 -3.24413
H 14.10338 1.40206 -4.13513
H 11.02971 5.89015 -0.82623
H 10.02898 6.89992 0.16274
H 9.45788 6.51817 -1.47887
H 7.94098 7.64744 0.40122
H 6.09572 7.74403 0.22124
H 7.14331 6.97649 -1.00629
H 4.78688 1.74183 3.20514
H 6.01379 -0.59200 4.58524
H 5.03022 2.53318 5.31314
H 6.80697 2.51729 5.03222
H 6.06674 1.17617 5.74573
H 3.59356 -0.47605 3.87754
H 4.41053 -1.97606 3.48068
H 3.54445 0.29746 1.84877
H 3.33044 -1.45498 1.49104
H 4.91544 -0.60794 1.20017
H 8.62390 -5.50896 1.61408
H 7.67210 -4.55337 2.79608
H 9.44911 -4.76213 3.06236
H 13.01885 -2.85957 -0.09071
H 14.62219 -4.61033 -3.68891
H 14.39930 -2.81478 -4.18918
H 13.15826 -4.08416 -4.33320

