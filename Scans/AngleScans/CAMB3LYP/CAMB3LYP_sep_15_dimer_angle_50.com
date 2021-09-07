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
Mg 7.88998 0.70806 0.60676
C 8.00017 -2.87836 0.67417
C 6.77320 0.58538 3.82737
C 7.35298 3.83047 0.43470
C 7.90347 0.50461 -2.87234
N 7.48028 -0.98711 1.96720
C 7.61754 -2.30297 1.90379
C 7.40629 -3.03753 3.25484
C 6.48719 -1.94078 3.98371
C 6.96235 -0.69246 3.24304
C 5.00321 -2.12203 3.71112
C 8.75232 -3.26474 4.02851
C 8.66257 -4.24345 5.20148
H 9.80828 -5.10572 5.47773
N 7.64344 2.05010 2.04784
C 7.20309 1.82432 3.30338
C 7.18024 3.13849 3.96413
C 7.48574 4.18054 2.94652
C 7.48541 3.37845 1.73957
C 6.96690 3.30317 5.44088
C 7.72440 5.65705 3.04212
O 7.73691 6.40585 2.00193
C 7.85002 6.21450 4.39673
N 7.44781 1.97739 -0.97932
C 7.31436 3.30166 -0.78849
C 7.27290 4.04098 -2.10457
C 7.65081 2.95723 -3.21824
C 7.76307 1.70341 -2.32341
C 5.93614 4.69120 -2.46050
C 8.83638 3.30098 -4.16107
C 10.12951 3.32033 -3.37823
N 8.12708 -0.86187 -0.85856
C 8.07348 -0.73234 -2.26371
C 8.17437 -2.09285 -2.83051
C 8.25717 -2.92813 -1.70626
C 8.19080 -2.14530 -0.52209
C 8.13285 -2.42261 -4.21458
C 8.37492 -4.29502 -1.25188
O 8.55667 -5.38059 -1.80163
C 8.08649 -4.33474 0.31069
C 9.16528 -5.15651 0.96141
O 10.34580 -4.76538 0.95926
O 8.59667 -6.15368 1.66333
C 9.62596 -6.99013 2.33619
H 6.20476 0.65251 4.75540
H 7.27649 4.91382 0.31647
H 8.01355 0.36525 -3.94990
H 6.87485 -3.97882 3.18188
H 6.72257 -1.86874 5.04238
H 4.81977 -2.88994 2.97437
H 4.65043 -1.16328 3.35422
H 4.46845 -2.31081 4.64639
H 9.01367 -2.34631 4.56972
H 9.52968 -3.57505 3.35765
H 7.73643 -4.83054 5.17851
H 8.59275 -3.66885 6.12829
H 6.60456 2.37553 5.90087
H 6.13994 4.01023 5.56761
H 7.84847 3.69499 5.95271
H 6.94556 5.99977 4.96091
H 8.10605 7.26657 4.25720
H 8.67083 5.70232 4.91108
H 8.00775 4.83640 -2.12462
H 6.81531 2.77568 -3.91020
H 6.01530 5.77931 -2.64442
H 5.23542 4.49069 -1.66092
H 5.42987 4.26328 -3.32408
H 8.71068 4.30342 -4.57757
H 8.93699 2.61563 -4.99603
H 10.13839 4.32331 -2.94114
H 10.98021 3.17063 -4.04000
H 10.10003 2.60083 -2.55174
H 8.90786 -3.15854 -4.38268
H 8.27751 -1.54238 -4.83752
H 7.11594 -2.81677 -4.29075
H 7.09758 -4.77310 0.45025
H 9.45752 -8.02300 2.01604
H 9.52605 -6.77505 3.41670
H 10.65943 -6.84417 2.04451

