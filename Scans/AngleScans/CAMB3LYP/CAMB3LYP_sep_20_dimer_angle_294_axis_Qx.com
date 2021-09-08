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
Mg 10.51697 0.94245 0.80621
C 8.04637 0.14966 3.28546
C 11.90480 3.06816 3.08434
C 12.43294 1.99661 -1.49275
C 8.32406 -0.32921 -1.58539
N 9.95733 1.45572 2.88362
C 9.03574 1.04536 3.74188
C 9.22981 1.56539 5.19152
C 10.09027 2.88065 4.86199
C 10.72261 2.43542 3.54495
C 9.23490 4.10704 4.59061
C 10.03764 0.55926 6.08434
C 9.99917 0.85701 7.58488
H 9.98513 -0.26674 8.51733
N 12.22860 1.94065 0.91794
C 12.64110 2.73840 1.92527
C 13.95762 3.25297 1.51764
C 14.21507 2.82485 0.11581
C 12.93181 2.24232 -0.22180
C 14.86505 3.99837 2.45272
C 15.41324 2.93452 -0.77782
O 15.33702 2.70810 -2.03708
C 16.65728 3.40530 -0.15154
N 10.31577 1.04476 -1.26052
C 11.32000 1.49796 -2.03129
C 11.05364 1.22239 -3.49194
C 9.76912 0.27020 -3.52223
C 9.44768 0.23076 -2.01209
C 10.79693 2.45077 -4.36462
C 9.90911 -1.07645 -4.28355
C 10.87766 -1.98176 -3.55710
N 8.62524 -0.10106 0.82912
C 7.86773 -0.51838 -0.28717
C 6.60342 -1.08071 0.23027
C 6.70402 -0.92932 1.62147
C 7.93472 -0.29758 1.94671
C 5.53608 -1.59335 -0.55974
C 6.04428 -1.15064 2.88814
O 5.01650 -1.72613 3.24342
C 6.82828 -0.34725 4.01344
C 7.03856 -1.27050 5.18216
O 7.76552 -2.27355 5.07243
O 6.52647 -0.70322 6.28955
C 6.72052 -1.57635 7.47774
H 12.28366 3.92685 3.63933
H 13.10501 2.24134 -2.31852
H 7.63088 -0.83522 -2.26076
H 8.31435 1.83836 5.70258
H 10.85602 3.05046 5.61442
H 8.18368 3.86362 4.54786
H 9.57216 4.49616 3.63882
H 9.44874 4.88638 5.32765
H 11.11027 0.69292 5.89371
H 9.71711 -0.44904 5.90768
H 9.21487 1.57938 7.84179
H 10.92846 1.36021 7.86275
H 14.33826 4.30601 3.36445
H 15.12791 4.94062 1.95985
H 15.78191 3.44908 2.67703
H 16.49277 4.38780 0.28445
H 17.42124 3.33204 -0.92799
H 16.90911 2.73213 0.67553
H 11.88792 0.70158 -3.94563
H 8.91315 0.75618 -4.01312
H 11.49293 2.53100 -5.22091
H 10.85910 3.33211 -3.74020
H 9.78981 2.52618 -4.77109
H 10.32882 -0.90730 -5.27823
H 8.96692 -1.59775 -4.41596
H 11.85280 -1.63665 -3.91361
H 10.70474 -3.02195 -3.82568
H 10.84214 -1.80882 -2.47521
H 5.20846 -2.50831 -0.08434
H 5.85001 -1.76832 -1.58668
H 4.82310 -0.76806 -0.48341
H 6.22183 0.51917 4.28012
H 5.73232 -1.74044 7.91845
H 7.47321 -1.06351 8.10560
H 7.04608 -2.59518 7.30323

