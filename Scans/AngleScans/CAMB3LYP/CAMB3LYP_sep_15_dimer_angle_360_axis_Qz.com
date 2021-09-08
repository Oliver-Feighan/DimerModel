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
Mg 7.88947 0.71461 0.60577
C 7.66863 -1.50811 3.41465
C 6.63278 3.15618 2.62915
C 7.64762 2.57345 -1.95424
C 8.02943 -2.10504 -1.43768
N 7.27022 0.70331 2.72913
C 7.28921 -0.18287 3.71334
C 6.95556 0.40200 5.11197
C 6.11090 1.68411 4.64121
C 6.72860 1.88992 3.25963
C 4.62924 1.39670 4.46396
C 8.24176 0.82169 5.90668
C 8.01412 1.11105 7.39201
H 9.06373 0.74892 8.34051
N 7.70667 2.68555 0.46242
C 7.19554 3.52478 1.38747
C 7.26531 4.86887 0.79357
C 7.70691 4.73422 -0.62120
C 7.68365 3.29275 -0.76855
C 7.00649 6.12067 1.58061
C 8.07527 5.73738 -1.67200
O 8.19952 5.40704 -2.90418
C 8.19468 7.13482 -1.23110
N 7.63166 0.30360 -1.41585
C 7.61191 1.29354 -2.32555
C 7.69310 0.74523 -3.73024
C 8.01655 -0.81227 -3.56545
C 7.97636 -0.91748 -2.02509
C 6.43735 0.91633 -4.58476
C 9.26686 -1.35413 -4.31096
C 10.52256 -0.76979 -3.70503
N 8.04306 -1.41920 0.90796
C 8.06018 -2.42237 -0.08574
C 8.05991 -3.72517 0.61078
C 8.01908 -3.38707 1.97185
C 7.97517 -1.97370 2.11311
C 8.04621 -5.00332 -0.01553
C 7.99241 -3.90465 3.32096
O 8.09699 -5.02218 3.82487
C 7.63663 -2.71404 4.31205
C 8.60774 -2.75921 5.45984
O 9.81811 -2.54348 5.27256
O 7.92162 -2.83348 6.61487
C 8.84122 -2.86955 7.78314
H 6.03463 3.93100 3.10952
H 7.67542 3.17030 -2.86883
H 8.17117 -3.02945 -2.00164
H 6.34386 -0.23705 5.73744
H 6.30751 2.54236 5.27848
H 4.40730 0.34521 4.57004
H 4.38072 1.73716 3.46728
H 4.04092 2.01413 5.14879
H 8.56312 1.81532 5.56879
H 9.01482 0.08631 5.79613
H 7.03998 0.74521 7.73865
H 7.95848 2.19357 7.53011
H 6.54204 5.89805 2.54918
H 6.24302 6.68789 1.03725
H 7.89805 6.74200 1.68895
H 7.25164 7.45837 -0.79699
H 8.55135 7.68683 -2.10285
H 8.94309 7.18697 -0.43249
H 8.49774 1.21471 -4.28273
H 7.19755 -1.44087 -3.94484
H 6.62319 1.46138 -5.52946
H 5.68850 1.42567 -3.99287
H 5.93057 -0.00954 -4.85148
H 9.25075 -1.03792 -5.35686
H 9.33917 -2.43658 -4.30286
H 10.60479 0.20338 -4.19825
H 11.38289 -1.39828 -3.92610
H 10.39302 -0.58571 -2.63219
H 8.75708 -5.61949 0.51879
H 8.29654 -4.93114 -1.07189
H 7.00161 -5.28574 0.14041
H 6.60685 -2.85836 4.64114
H 8.59263 -3.76712 8.35782
H 8.71640 -1.89515 8.29178
H 9.89497 -3.02925 7.58678

