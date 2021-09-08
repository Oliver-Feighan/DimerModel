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
Mg 8.94153 0.80925 0.67278
C 8.35671 4.23079 1.58387
C 8.22244 1.69814 -2.54085
C 9.00765 -2.27021 -0.08879
C 8.49179 -0.01611 4.02870
N 8.37730 2.74524 -0.23596
C 8.24624 3.98957 0.19857
C 8.06155 5.04649 -0.92320
C 7.46988 4.08453 -2.06489
C 8.08595 2.76188 -1.61352
C 5.95716 3.94907 -2.01331
C 9.42122 5.69323 -1.36482
C 9.28545 6.94695 -2.23172
H 10.26567 8.01673 -2.06677
N 9.13890 -0.07281 -1.09403
C 8.81898 0.44092 -2.30026
C 9.13370 -0.60931 -3.28114
C 9.51069 -1.84434 -2.54135
C 9.20697 -1.44101 -1.18298
C 9.13784 -0.36461 -4.76220
C 10.04287 -3.17345 -2.98456
O 10.07474 -4.18307 -2.19558
C 10.43871 -3.28347 -4.39619
N 8.56615 -0.92107 1.76248
C 8.71805 -2.13542 1.20553
C 8.66198 -3.22518 2.24941
C 8.68271 -2.47090 3.65945
C 8.65598 -1.00817 3.16440
C 7.44473 -4.14734 2.18368
C 9.78894 -2.89053 4.66586
C 11.14637 -2.48189 4.14072
N 8.68752 1.89784 2.52169
C 8.48975 1.35567 3.81050
C 8.25336 2.48909 4.72806
C 8.30758 3.61952 3.89880
C 8.53987 3.21750 2.55570
C 7.98023 2.38843 6.12139
C 8.21057 5.06075 3.85253
O 8.10861 5.95202 4.69454
C 8.11216 5.51205 2.33166
C 9.08071 6.64310 2.11918
O 10.30568 6.45362 2.22000
O 8.41802 7.70309 1.62154
C 9.33823 8.84692 1.38346
H 7.79522 1.81989 -3.53661
H 9.13035 -3.33904 -0.27863
H 8.44065 -0.18350 5.10668
H 7.35171 5.83236 -0.69449
H 7.84178 4.36269 -3.04758
H 5.53910 4.42994 -1.14146
H 5.75763 2.88570 -1.99110
H 5.51363 4.31961 -2.94192
H 9.92071 5.02488 -2.07796
H 10.03555 5.90826 -0.51218
H 8.26672 7.35287 -2.21273
H 9.44241 6.66553 -3.27584
H 8.66003 0.59036 -5.01346
H 8.48652 -1.12114 -5.21311
H 10.13452 -0.44980 -5.20041
H 9.58520 -3.05365 -5.02956
H 10.87632 -4.27823 -4.50030
H 11.20007 -2.52398 -4.60625
H 9.52999 -3.86955 2.18277
H 7.75001 -2.62992 4.22036
H 7.71185 -5.21613 2.08165
H 6.82090 -3.83083 1.35833
H 6.76346 -4.07383 3.02965
H 9.81167 -3.97790 4.77199
H 9.65147 -2.47296 5.65769
H 11.40400 -3.29910 3.46036
H 11.86455 -2.40458 4.95453
H 11.07772 -1.56513 3.54367
H 8.57033 3.15214 6.61030
H 8.21711 1.39721 6.50238
H 6.90437 2.58281 6.12379
H 7.08115 5.81563 2.14608
H 8.93370 9.70223 1.93345
H 9.41441 8.94527 0.28425
H 10.33701 8.78154 1.79893

