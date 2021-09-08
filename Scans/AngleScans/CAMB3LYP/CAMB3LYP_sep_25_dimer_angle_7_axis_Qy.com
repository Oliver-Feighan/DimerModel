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
Mg 13.14578 1.18772 1.00409
C 13.24970 -1.00517 3.84302
C 12.18540 3.66614 3.14177
C 12.60927 3.02584 -1.52591
C 12.99083 -1.65267 -1.00934
N 12.79777 1.20530 3.18831
C 12.92836 0.32819 4.17217
C 12.78064 0.93091 5.59500
C 11.90065 2.22071 5.21961
C 12.34278 2.40457 3.76930
C 10.40478 1.95323 5.23267
C 14.16184 1.33927 6.21752
C 14.12624 1.64580 7.71636
H 15.28168 1.27733 8.52979
N 12.97258 3.15968 0.86289
C 12.59274 4.01488 1.83533
C 12.60524 5.35219 1.22248
C 12.86400 5.19790 -0.23490
C 12.80333 3.75560 -0.36215
C 12.46385 6.61493 2.02179
C 13.11086 6.18571 -1.33466
O 13.07512 5.84208 -2.56895
C 13.30317 7.58533 -0.92776
N 12.63091 0.76160 -0.96458
C 12.51029 1.74316 -1.87557
C 12.40730 1.18061 -3.27318
C 12.72819 -0.37982 -3.13294
C 12.88022 -0.47002 -1.59862
C 11.05663 1.36195 -3.96542
C 13.86778 -0.94677 -4.02322
C 15.19722 -0.37511 -3.58598
N 13.30778 -0.94518 1.30833
C 13.18677 -1.95775 0.33153
C 13.25663 -3.25383 1.03702
C 13.39140 -2.90245 2.38863
C 13.38430 -1.48733 2.51857
C 13.14751 -4.53745 0.43161
C 13.52737 -3.40694 3.73611
O 13.67953 -4.52111 4.23533
C 13.31458 -2.20206 4.75069
C 14.42136 -2.25060 5.76814
O 15.60141 -2.05426 5.42825
O 13.88467 -2.30407 7.00082
C 14.94304 -2.34258 8.04494
H 11.66258 4.45404 3.68469
H 12.53000 3.61363 -2.44335
H 13.04842 -2.58428 -1.57631
H 12.24384 0.30671 6.29927
H 12.18703 3.08194 5.81764
H 10.18399 0.90612 5.37741
H 10.03770 2.28793 4.27128
H 9.91528 2.58554 5.97894
H 14.45140 2.32491 5.83100
H 14.90510 0.59172 6.01913
H 13.19853 1.29743 8.18641
H 14.10272 2.73026 7.84831
H 12.12169 6.40817 3.04335
H 11.64582 7.18809 1.57214
H 13.37011 7.22422 2.01062
H 12.42642 7.92662 -0.38250
H 13.55491 8.12392 -1.84342
H 14.14650 7.63407 -0.22988
H 13.14240 1.63314 -3.92735
H 11.85980 -0.99997 -3.39968
H 11.12966 1.89537 -4.93196
H 10.39478 1.88764 -3.29003
H 10.50815 0.44109 -4.15621
H 13.72475 -0.64016 -5.06230
H 13.92618 -2.03003 -4.01222
H 15.22981 0.59211 -4.09640
H 16.01460 -1.01808 -3.90615
H 15.20578 -0.17913 -2.50748
H 13.91158 -5.15887 0.87943
H 13.26424 -4.47881 -0.64854
H 12.12708 -4.80318 0.72037
H 12.33241 -2.32830 5.20784
H 14.75665 -3.23102 8.65618
H 14.89596 -1.36175 8.55436
H 15.96161 -2.51942 7.71984

