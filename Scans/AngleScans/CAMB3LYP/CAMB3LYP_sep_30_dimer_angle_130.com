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
Mg 15.77424 1.41241 1.19728
C 15.89658 0.71666 -2.32124
C 15.15317 -1.90658 1.68021
C 15.21463 2.06326 4.25190
C 15.26601 4.76106 0.37631
N 15.57475 -0.27665 -0.21677
C 15.70215 -0.43190 -1.52597
C 15.69657 -1.90936 -2.00202
C 14.89616 -2.55620 -0.76921
C 15.25372 -1.54074 0.31414
C 13.38826 -2.54293 -0.95733
C 17.14362 -2.50392 -2.12271
C 17.23169 -3.83427 -2.87383
H 18.40665 -4.08915 -3.70266
N 15.74536 0.21650 2.78080
C 15.49848 -1.11021 2.79416
C 15.57378 -1.51926 4.20525
C 15.72225 -0.29671 5.04075
C 15.54166 0.73238 4.03654
C 15.58416 -2.95718 4.63618
C 15.97119 -0.08828 6.50389
O 15.82688 1.06009 7.05457
C 16.29803 -1.28657 7.29051
N 15.09804 3.11833 2.17450
C 14.99350 3.15175 3.51460
C 14.75446 4.55875 4.00817
C 14.96209 5.50222 2.73367
C 15.20842 4.42259 1.65705
C 13.37885 4.82467 4.61920
C 15.99254 6.65613 2.87265
C 17.38843 6.08929 2.99750
N 15.79034 2.58942 -0.61447
C 15.52651 3.97131 -0.73627
C 15.54251 4.29217 -2.17822
C 15.79374 3.06054 -2.80147
C 15.90498 2.04027 -1.81853
C 15.29420 5.57257 -2.74834
C 15.97956 2.38966 -4.06800
O 16.07780 2.75545 -5.23860
C 15.92876 0.82036 -3.82085
C 17.09369 0.19873 -4.54130
O 18.26017 0.44565 -4.18781
O 16.63768 -0.74758 -5.38209
C 17.75701 -1.40069 -6.11174
H 14.73027 -2.88091 1.92661
H 15.12028 2.36222 5.29841
H 15.21338 5.80009 0.04417
H 15.16107 -2.08656 -2.92702
H 15.28758 -3.53798 -0.51569
H 13.09710 -1.99197 -1.83926
H 12.98506 -2.07513 -0.06871
H 12.99997 -3.56545 -0.96223
H 17.48233 -2.82609 -1.12962
H 17.81186 -1.79267 -2.56790
H 16.31311 -4.05452 -3.43120
H 17.30111 -4.64230 -2.14153
H 15.29574 -3.62522 3.81525
H 14.78490 -3.07508 5.37588
H 16.53215 -3.25411 5.08984
H 15.48860 -2.00749 7.20357
H 16.52927 -0.92344 8.29378
H 17.18716 -1.75788 6.85700
H 15.47726 4.83014 4.76778
H 14.03243 6.01902 2.45352
H 13.42853 5.21137 5.65461
H 12.80614 3.90749 4.58289
H 12.74919 5.51259 4.05750
H 15.80488 7.22329 3.78767
H 15.96742 7.35942 2.04697
H 17.46183 5.84633 4.06177
H 18.13040 6.83148 2.71000
H 17.48389 5.15150 2.43796
H 16.03592 5.71784 -3.52246
H 15.34300 6.35839 -1.99742
H 14.27773 5.42542 -3.12305
H 14.97237 0.45959 -4.20129
H 17.54342 -1.30040 -7.18038
H 17.82008 -2.42744 -5.70485
H 18.73492 -0.93811 -6.04719

