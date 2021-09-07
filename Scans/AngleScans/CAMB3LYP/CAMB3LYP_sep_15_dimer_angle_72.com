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
Mg 7.89032 0.70549 0.60539
C 8.08074 -2.64714 -0.66058
C 6.90160 -0.63409 3.58240
C 7.27951 3.65098 1.61458
C 7.77392 1.79946 -2.70148
N 7.57028 -1.38481 1.25384
C 7.73347 -2.57907 0.70478
C 7.59050 -3.76674 1.69400
C 6.67666 -3.04822 2.80192
C 7.09570 -1.59923 2.56213
C 5.18755 -3.16716 2.52311
C 8.97001 -4.21670 2.29109
C 8.94683 -5.56133 3.02129
H 10.12077 -6.42475 2.92725
N 7.67055 1.41246 2.44671
C 7.28398 0.72475 3.54165
C 7.25817 1.70118 4.64169
C 7.50150 3.05481 4.07318
C 7.47197 2.75462 2.65575
C 7.09847 1.30245 6.08007
C 7.71161 4.39938 4.70103
O 7.66773 5.47883 4.01133
C 7.87731 4.42218 6.16168
N 7.35985 2.45408 -0.38626
C 7.20518 3.60932 0.28414
C 7.09694 4.77982 -0.66380
C 7.45499 4.19661 -2.10923
C 7.62885 2.70589 -1.74466
C 5.73366 5.46905 -0.71681
C 8.59556 4.90427 -2.89089
C 9.91721 4.67815 -2.19269
N 8.10465 -0.20457 -1.34240
C 7.99406 0.43189 -2.59810
C 8.10247 -0.61953 -3.63017
C 8.24669 -1.80705 -2.89702
C 8.20912 -1.51864 -1.50601
C 8.01472 -0.41708 -5.03639
C 8.41150 -3.24025 -2.98347
O 8.59538 -4.03977 -3.90036
C 8.18453 -3.86307 -1.53884
C 9.30519 -4.82921 -1.26825
O 10.47598 -4.42447 -1.15858
O 8.78584 -6.03384 -0.96903
C 9.85822 -7.02347 -0.68192
H 6.36807 -0.93339 4.48498
H 7.17505 4.69830 1.90730
H 7.84532 2.07100 -3.75686
H 7.07718 -4.63253 1.29315
H 6.95111 -3.36345 3.80525
H 4.99249 -3.61522 1.56012
H 4.80055 -2.15710 2.55580
H 4.69350 -3.70568 3.33689
H 9.23208 -3.55404 3.12592
H 9.72743 -4.23091 1.53166
H 8.03344 -6.13009 2.80878
H 8.90039 -5.37161 4.09640
H 6.77434 0.25869 6.17454
H 6.26188 1.88411 6.48213
H 7.99044 1.50808 6.67570
H 7.00017 3.98363 6.63153
H 8.10489 5.45974 6.41381
H 8.72827 3.78503 6.42723
H 7.81305 5.55139 -0.40892
H 6.59754 4.25429 -2.79572
H 5.78205 6.55036 -0.48769
H 5.06883 4.96394 -0.02876
H 5.20383 5.37247 -1.66294
H 8.43217 5.98468 -2.90372
H 8.67869 4.57886 -3.92249
H 9.92121 5.44904 -1.41628
H 10.74478 4.81230 -2.88635
H 9.93526 3.70409 -1.69017
H 8.79844 -1.01209 -5.48611
H 8.11612 0.63520 -5.29366
H 7.00442 -0.79014 -5.22447
H 7.21148 -4.35574 -1.54369
H 9.69999 -7.87072 -1.35632
H 9.79546 -7.22542 0.40391
H 10.87626 -6.74481 -0.92768

