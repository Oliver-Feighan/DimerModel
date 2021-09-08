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
Mg 7.88503 0.71450 0.60809
C 6.48322 -1.66151 2.90342
C 5.65851 2.98866 1.83506
C 8.78565 2.69427 -1.70207
C 9.16842 -1.98411 -1.18526
N 6.32375 0.56383 2.16757
C 5.92368 -0.37164 3.01561
C 4.92262 0.12324 4.09382
C 4.32863 1.39075 3.30675
C 5.51975 1.69608 2.40083
C 3.13609 1.04928 2.42893
C 5.63909 0.55548 5.42116
C 4.70786 0.75554 6.61888
H 5.19313 0.38939 7.94663
N 7.66889 2.68051 0.44478
C 6.72401 3.44659 1.02925
C 6.98509 4.82205 0.57744
C 8.05975 4.78185 -0.45126
C 8.20173 3.35057 -0.62836
C 6.30103 6.01797 1.17346
C 8.82320 5.85432 -1.16760
O 9.54495 5.59548 -2.19468
C 8.62705 7.23027 -0.68825
N 8.65709 0.40051 -1.29656
C 9.01410 1.43551 -2.07709
C 9.79493 0.96649 -3.28155
C 10.09755 -0.58210 -3.02099
C 9.32887 -0.77048 -1.69470
C 9.09633 1.12815 -4.63143
C 11.58381 -1.02903 -3.08383
C 12.35386 -0.42384 -1.93224
N 8.00963 -1.42368 0.89217
C 8.56581 -2.37106 0.00492
C 8.31358 -3.70775 0.58117
C 7.60236 -3.44439 1.76155
C 7.40634 -2.04372 1.90016
C 8.68372 -4.95028 -0.00640
C 6.96360 -4.03345 2.91640
O 6.88388 -5.17055 3.37941
C 6.10056 -2.91359 3.64267
C 6.40128 -2.97749 5.11509
O 7.53626 -2.69994 5.54100
O 5.25083 -3.14253 5.79302
C 5.49608 -3.20067 7.25873
H 4.85530 3.71029 1.98661
H 9.21156 3.33943 -2.47398
H 9.62209 -2.87031 -1.63415
H 4.12760 -0.57395 4.32967
H 4.13988 2.22170 3.98147
H 2.95775 -0.01497 2.38774
H 3.37604 1.43104 1.44519
H 2.25321 1.60348 2.76018
H 6.01946 1.57860 5.30595
H 6.41508 -0.13885 5.67879
H 3.71246 0.33001 6.44267
H 4.52410 1.82579 6.74061
H 5.44349 5.72441 1.79140
H 5.85834 6.57963 0.34375
H 6.98935 6.67069 1.71462
H 7.57305 7.48923 -0.75524
H 9.32288 7.84266 -1.26504
H 8.89459 7.27235 0.37348
H 10.73441 1.49899 -3.36471
H 9.60341 -1.22455 -3.76457
H 9.67823 1.73006 -5.35473
H 8.12456 1.57249 -4.46181
H 8.84001 0.19666 -5.13326
H 12.05223 -0.65880 -3.99907
H 11.71189 -2.10624 -3.06952
H 12.60101 0.57676 -2.29940
H 13.25251 -1.00190 -1.72650
H 11.71333 -0.30270 -1.05094
H 9.08786 -5.56257 0.78876
H 9.40571 -4.81142 -0.80840
H 7.71303 -5.28552 -0.38143
H 5.05086 -3.11959 3.42990
H 5.05948 -4.13736 7.61881
H 5.08062 -2.26093 7.66865
H 6.52226 -3.30397 7.59148

