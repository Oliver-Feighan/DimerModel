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
Mg 31.54003 2.84243 2.38874
C 30.99552 6.38140 2.63011
C 30.88527 3.12010 -0.94724
C 31.56908 -0.32589 2.22048
C 31.02505 2.66873 5.83118
N 31.02083 4.57996 1.12238
C 30.90280 5.88537 1.31299
C 30.75241 6.71470 0.00948
C 30.16338 5.56277 -0.94182
C 30.75127 4.34091 -0.23887
C 28.64808 5.45860 -0.89337
C 32.12889 7.24956 -0.52083
C 32.02656 8.31948 -1.61012
H 33.02096 9.38869 -1.63116
N 31.75085 1.64140 0.82297
C 31.45794 1.92311 -0.46393
C 31.77117 0.70324 -1.22417
C 32.11694 -0.37529 -0.25880
C 31.79855 0.28013 0.99385
C 31.80225 0.66495 -2.72450
C 32.63478 -1.77068 -0.43468
O 32.63832 -2.61421 0.53034
C 33.05075 -2.14923 -1.79305
N 31.12029 1.35284 3.77701
C 31.26154 0.05358 3.46094
C 31.17193 -0.81960 4.68972
C 31.18269 0.18608 5.93313
C 31.18690 1.52982 5.17174
C 29.94136 -1.72217 4.77614
C 32.26631 -0.05076 7.02049
C 33.63807 0.23458 6.45285
N 31.27459 4.26246 3.99532
C 31.04820 3.97490 5.35919
C 30.81558 5.26357 6.04305
C 30.90066 6.21709 5.01736
C 31.14740 5.56673 3.77815
C 30.51924 5.43022 7.42527
C 30.82727 7.62501 4.69944
O 30.72639 8.65995 5.35705
C 30.75973 7.78343 3.11933
C 31.74931 8.84197 2.71593
O 32.96940 8.65932 2.87291
O 31.11137 9.79776 2.01602
C 32.05322 10.86465 1.58420
H 30.47558 3.05781 -1.95578
H 31.67773 -1.41281 2.23705
H 30.95449 2.70773 6.92029
H 30.05166 7.53849 0.07350
H 30.55490 5.64641 -1.95230
H 28.22418 6.10012 -0.13514
H 28.43136 4.42101 -0.67548
H 28.22500 5.65347 -1.88299
H 32.62875 6.45274 -1.08649
H 32.73321 7.61331 0.28732
H 31.01423 8.73461 -1.68634
H 32.19526 7.84475 -2.57976
H 31.34364 1.56161 -3.15935
H 31.14610 -0.15456 -3.03711
H 32.80416 0.48624 -3.12062
H 32.21096 -2.03183 -2.47380
H 33.47407 -3.15128 -1.70044
H 33.82725 -1.45251 -2.12808
H 32.03053 -1.47596 4.76118
H 30.23896 0.14722 6.49678
H 30.19303 -2.79437 4.88158
H 29.33556 -1.55866 3.89477
H 29.24825 -1.48224 5.58065
H 32.27012 -1.09897 7.32937
H 32.12007 0.54760 7.91356
H 33.89325 -0.69920 5.94291
H 34.34463 0.45447 7.25063
H 33.59329 1.02349 5.69306
H 31.11370 6.26472 7.77273
H 30.73438 4.52543 7.98992
H 29.44669 5.63520 7.37145
H 29.73669 8.05973 2.86122
H 31.65382 11.81319 1.95625
H 32.14806 10.75353 0.48765
H 33.04425 10.86592 2.02273

