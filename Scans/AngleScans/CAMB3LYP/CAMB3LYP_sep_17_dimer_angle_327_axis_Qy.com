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
Mg 8.93550 0.80909 0.68806
C 7.40924 -1.58864 2.87905
C 6.62016 3.06145 1.78382
C 9.96531 2.80727 -1.55119
C 10.34819 -1.87110 -1.03431
N 7.28046 0.63946 2.14557
C 6.83458 -0.30323 2.96219
C 5.76427 0.17907 3.97776
C 5.21188 1.44669 3.16115
C 6.45535 1.76496 2.33339
C 4.07929 1.10182 2.20845
C 6.39286 0.60906 5.34962
C 5.38683 0.79638 6.48741
H 5.79033 0.42644 7.84122
N 8.71604 2.77443 0.52114
C 7.73088 3.53093 1.04882
C 8.01002 4.91052 0.62107
C 9.14742 4.88314 -0.33824
C 9.31044 3.45382 -0.51308
C 7.28136 6.09856 1.17879
C 9.94670 5.96460 -0.99991
O 10.73337 5.71616 -1.98085
C 9.71103 7.33660 -0.52708
N 9.82788 0.51054 -1.16582
C 10.22584 1.55211 -1.91728
C 11.08411 1.09486 -3.07257
C 11.38078 -0.45300 -2.80111
C 10.53168 -0.65370 -1.52669
C 10.47053 1.25898 -4.46289
C 12.87120 -0.88948 -2.77261
C 13.56307 -0.28523 -1.57195
N 9.05724 -1.32968 0.96891
C 9.67479 -2.26850 0.11373
C 9.39639 -3.60994 0.66642
C 8.61058 -3.35773 1.80104
C 8.39626 -1.95920 1.93394
C 9.81155 -4.84676 0.09718
C 7.90475 -3.95729 2.91056
O 7.80421 -5.09737 3.36205
C 6.98984 -2.84723 3.58664
C 7.19792 -2.91695 5.07471
O 8.30190 -2.63397 5.57248
O 6.00837 -3.09344 5.67819
C 6.16147 -3.15774 7.15611
H 5.80389 3.77679 1.88813
H 10.43427 3.45942 -2.29163
H 10.83547 -2.75178 -1.45815
H 4.96099 -0.52477 4.15976
H 4.97520 2.27272 3.82675
H 3.91148 0.03661 2.15091
H 4.37783 1.49045 1.24363
H 3.17343 1.64822 2.48627
H 6.77241 1.63535 5.26357
H 7.15606 -0.08134 5.65210
H 4.40752 0.36503 6.24690
H 5.18816 1.86468 6.60260
H 6.38881 5.79587 1.74016
H 6.88767 6.66163 0.32566
H 7.92966 6.75305 1.76533
H 8.66150 7.58873 -0.65891
H 10.43733 7.95678 -1.05599
H 9.91105 7.37481 0.54955
H 12.02314 1.63418 -3.09390
H 10.93892 -1.09482 -3.57742
H 11.09242 1.86869 -5.14522
H 9.48690 1.69579 -4.35249
H 10.25289 0.32846 -4.98440
H 13.39354 -0.51118 -3.65477
H 13.00580 -1.96585 -2.75557
H 13.82565 0.71898 -1.91793
H 14.45111 -0.85825 -1.31297
H 12.86760 -0.17315 -0.73206
H 10.16930 -5.46053 0.91315
H 10.58148 -4.69870 -0.65717
H 8.86875 -5.18659 -0.33975
H 5.95707 -3.05923 3.30731
H 5.70980 -4.09930 7.48343
H 5.71439 -2.22305 7.54371
H 7.16544 -3.25584 7.55218

