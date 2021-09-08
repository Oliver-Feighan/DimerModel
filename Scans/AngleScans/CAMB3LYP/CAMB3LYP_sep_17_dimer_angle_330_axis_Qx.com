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
Mg 8.94173 0.80536 0.68609
C 7.51779 -1.11371 3.36352
C 8.90115 3.49081 2.78878
C 9.77155 2.62245 -1.77907
C 7.82468 -1.67070 -1.49726
N 8.27273 1.03359 2.78195
C 7.80644 0.22066 3.71792
C 7.71573 0.85340 5.13248
C 7.62082 2.39699 4.69981
C 8.33948 2.31836 3.35444
C 6.19663 2.86657 4.45383
C 8.99715 0.57624 5.99454
C 8.85200 0.89668 7.48384
H 9.54145 0.04530 8.44940
N 9.73843 2.62244 0.64068
C 9.64205 3.57715 1.58964
C 10.38392 4.73664 1.07063
C 10.78662 4.44615 -0.33228
C 10.08112 3.20031 -0.55645
C 10.71469 5.93441 1.91287
C 11.65186 5.17619 -1.31443
O 11.67231 4.86256 -2.55703
C 12.40342 6.32858 -0.79615
N 8.63448 0.62850 -1.36248
C 9.14576 1.53045 -2.21864
C 9.03357 1.05186 -3.64649
C 8.55760 -0.47203 -3.55269
C 8.38341 -0.58901 -2.02254
C 8.06649 1.83299 -4.53587
C 9.43336 -1.52889 -4.27986
C 10.77816 -1.64181 -3.59859
N 8.03198 -1.14468 0.88033
C 7.62129 -2.00215 -0.16378
C 6.95419 -3.16205 0.46235
C 7.00304 -2.88565 1.83702
C 7.63661 -1.63177 2.05118
C 6.36313 -4.25581 -0.23080
C 6.65322 -3.36449 3.15499
O 6.17805 -4.40740 3.60253
C 6.85802 -2.17929 4.19410
C 7.61987 -2.72112 5.37245
O 8.79319 -3.11184 5.24100
O 6.91766 -2.48781 6.49621
C 7.63766 -2.99737 7.69361
H 8.72318 4.44401 3.28741
H 10.13556 3.15760 -2.65918
H 7.53623 -2.53184 -2.10381
H 6.83725 0.57180 5.70059
H 8.16920 3.03462 5.38833
H 5.49044 2.05069 4.49597
H 6.20027 3.31327 3.46818
H 5.93979 3.67142 5.14852
H 9.77577 1.30012 5.72162
H 9.32587 -0.43762 5.87343
H 7.80395 1.03731 7.77497
H 9.31630 1.86682 7.67691
H 10.14570 5.93616 2.85077
H 10.35098 6.81531 1.37288
H 11.78733 6.04414 2.08629
H 11.70928 7.05490 -0.38012
H 13.03107 6.66435 -1.62392
H 13.03734 5.98958 0.03079
H 9.99503 1.08975 -4.14373
H 7.56050 -0.61538 -3.99435
H 8.54554 2.24731 -5.44302
H 7.62255 2.62350 -3.94549
H 7.19299 1.27540 -4.86947
H 9.63145 -1.21431 -5.30746
H 8.97529 -2.51123 -4.32643
H 11.34667 -0.81580 -4.03630
H 11.24100 -2.60112 -3.82116
H 10.69184 -1.44902 -2.52293
H 6.65780 -5.15386 0.29567
H 6.67747 -4.28321 -1.27205
H 5.30453 -4.00264 -0.12812
H 5.86890 -1.81743 4.47741
H 6.95533 -3.67918 8.21056
H 7.96814 -2.09907 8.24822
H 8.49385 -3.64066 7.52753

