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
Mg 31.52278 2.84127 2.39207
C 30.78265 0.03138 0.28592
C 31.08842 4.82164 -0.35094
C 32.60640 5.28196 4.10560
C 32.98918 0.60396 4.62585
N 31.00337 2.37024 0.29430
C 30.68455 1.27313 -0.37580
C 30.15808 1.52786 -1.81368
C 30.84122 2.95843 -2.06997
C 30.95678 3.43666 -0.62407
C 32.24072 2.86104 -2.65414
C 28.59273 1.62059 -1.86939
C 27.99454 1.55910 -3.27656
H 26.71403 0.88140 -3.45952
N 31.34352 4.79594 2.09931
C 31.14816 5.42802 0.92306
C 31.10575 6.86393 1.24013
C 31.45647 7.03782 2.67591
C 31.84450 5.68420 3.01829
C 30.68026 7.90465 0.24550
C 31.45643 8.22831 3.58638
O 32.02584 8.19810 4.73429
C 30.85157 9.45902 3.05630
N 32.82060 2.96082 4.01185
C 33.08209 4.14003 4.60263
C 33.81758 3.94540 5.90707
C 33.78247 2.37166 6.18941
C 33.07960 1.89851 4.89810
C 35.26539 4.43481 5.93388
C 33.20318 1.91463 7.55624
C 31.72280 2.21512 7.61581
N 31.68195 0.69335 2.56290
C 32.36262 -0.03464 3.56331
C 32.28724 -1.46211 3.19037
C 31.58222 -1.46352 1.97734
C 31.25977 -0.12934 1.60928
C 32.86950 -2.54224 3.91176
C 31.04671 -2.28946 0.91916
O 30.93946 -3.50082 0.73250
C 30.61532 -1.34857 -0.28709
C 29.23803 -1.76156 -0.72871
O 28.26219 -1.61949 0.02903
O 29.26490 -2.05127 -2.04238
C 27.91923 -2.45571 -2.52964
H 31.19883 5.50792 -1.19100
H 32.91074 6.07845 4.78858
H 33.33668 -0.16421 5.32000
H 30.49865 0.81502 -2.55507
H 30.18501 3.61360 -2.63699
H 32.59199 1.84054 -2.69053
H 32.87181 3.45553 -2.00657
H 32.27444 3.34465 -3.63460
H 28.28338 2.63252 -1.57778
H 28.14419 0.86611 -1.25278
H 28.72290 1.21146 -4.01925
H 27.75146 2.57572 -3.59502
H 30.64188 7.49665 -0.77208
H 31.47717 8.65520 0.20733
H 29.74529 8.39574 0.52360
H 31.36908 9.75276 2.14616
H 30.86584 10.17139 3.88346
H 29.81461 9.24564 2.77381
H 33.31268 4.45891 6.71607
H 34.79102 1.93289 6.17709
H 35.46150 5.17403 6.73336
H 35.50283 4.85155 4.96409
H 36.01537 3.65274 6.03917
H 33.66691 2.47622 8.37093
H 33.36016 0.86092 7.76128
H 31.69639 3.26724 7.91487
H 31.23292 1.58271 8.35345
H 31.26634 2.14343 6.62175
H 32.13157 -3.33299 3.93771
H 33.16309 -2.23832 4.91436
H 33.73424 -2.75638 3.27800
H 31.35392 -1.47293 -1.07991
H 28.03072 -3.43932 -2.99630
H 27.57440 -1.62736 -3.17669
H 27.15738 -2.65942 -1.78630

