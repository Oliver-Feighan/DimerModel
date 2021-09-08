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
Mg 8.94075 0.80607 0.68685
C 8.90148 -2.40879 2.28128
C 7.72007 2.09860 3.59785
C 8.53797 3.56703 -0.82388
C 9.05425 -0.86542 -2.36912
N 8.41963 -0.13451 2.61972
C 8.50467 -1.35404 3.12950
C 8.22146 -1.43526 4.65356
C 7.32588 -0.11065 4.80265
C 7.87472 0.68950 3.62316
C 5.84493 -0.35566 4.56581
C 9.53228 -1.34190 5.51091
C 9.36593 -1.72277 6.98367
H 10.46611 -2.41119 7.65311
N 8.70472 2.64080 1.40560
C 8.21657 2.98399 2.61606
C 8.22711 4.45463 2.65191
C 8.60667 4.95403 1.30237
C 8.61102 3.71321 0.55374
C 7.97480 5.23983 3.90619
C 8.90264 6.32354 0.77037
O 8.97827 6.55526 -0.48797
C 9.00877 7.40341 1.76240
N 8.60067 1.28565 -1.30732
C 8.51581 2.56715 -1.70540
C 8.54575 2.67364 -3.21144
C 8.91338 1.20917 -3.73846
C 8.94614 0.45595 -2.39050
C 7.24854 3.14091 -3.87121
C 10.14085 1.08853 -4.68261
C 11.40867 1.40990 -3.92448
N 9.15874 -1.24475 0.04288
C 9.15430 -1.72711 -1.28429
C 9.21688 -3.20131 -1.21194
C 9.23026 -3.47752 0.16357
C 9.15921 -2.26226 0.89678
C 9.20500 -4.09002 -2.32377
C 9.27760 -4.52133 1.16203
O 9.43167 -5.74171 1.13630
C 8.93924 -3.88237 2.57753
C 9.96259 -4.37250 3.56497
O 11.15766 -4.04805 3.45028
O 9.33194 -4.96004 4.59814
C 10.30456 -5.45280 5.60958
H 7.12623 2.56952 4.38172
H 8.50973 4.49733 -1.39591
H 9.19204 -1.45464 -3.27825
H 7.65437 -2.30445 4.96452
H 7.53093 0.40142 5.73919
H 5.65317 -1.36021 4.21905
H 5.54309 0.36659 3.81861
H 5.27399 -0.11370 5.46681
H 9.81411 -0.28675 5.62040
H 10.31676 -1.92762 5.07257
H 8.41764 -2.24110 7.17080
H 9.29091 -0.80564 7.57295
H 7.56053 4.60679 4.70057
H 7.17401 5.95267 3.68154
H 8.85534 5.79175 4.24191
H 8.07912 7.47208 2.32230
H 9.31199 8.28843 1.19966
H 9.79144 7.14089 2.48290
H 9.31289 3.36640 -3.53497
H 8.09309 0.76917 -4.32444
H 7.37799 4.04372 -4.49751
H 6.51561 3.31821 -3.09538
H 6.75229 2.39707 -4.49221
H 10.06944 1.81938 -5.49188
H 10.23921 0.10983 -5.14029
H 11.44508 2.50294 -3.95663
H 12.27271 0.97146 -4.41970
H 11.32395 1.11370 -2.87246
H 9.95401 -4.84540 -2.12659
H 9.40500 -3.56430 -3.25511
H 8.17564 -4.45458 -2.27108
H 7.92930 -4.19531 2.84524
H 10.10393 -6.51906 5.75282
H 10.17996 -4.79413 6.48957
H 11.35174 -5.47018 5.33107

