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
Mg 7.87409 0.71203 0.59383
C 7.60424 -2.07485 -1.65106
C 7.91734 2.73070 -2.15530
C 8.56302 3.14327 2.51265
C 8.94546 -1.53475 3.03306
N 7.76744 0.26734 -1.57023
C 7.60426 -0.82399 -2.30283
C 7.35092 -0.55470 -3.81043
C 8.03801 0.89302 -3.91488
C 7.87001 1.34836 -2.46673
C 9.52398 0.83281 -4.22730
C 5.82205 -0.49131 -4.15763
C 5.49957 -0.54004 -5.65275
H 4.29137 -1.23911 -6.08176
N 7.70945 2.66763 0.29873
C 7.72388 3.31604 -0.88476
C 7.59092 4.74517 -0.56216
C 7.66252 4.90099 0.91612
C 8.00954 3.54942 1.30714
C 7.33625 5.79449 -1.60500
C 7.46547 6.07536 1.82620
O 7.81029 6.03635 3.06006
C 6.94341 7.30311 1.20850
N 8.84254 0.82867 2.42981
C 8.96245 2.00234 3.07484
C 9.44471 1.79947 4.49140
C 9.39232 0.22071 4.74113
C 8.95449 -0.24358 3.33467
C 10.85062 2.31625 4.79593
C 8.57756 -0.27097 5.96876
C 7.10591 -0.00023 5.75350
N 8.04615 -1.43503 0.76284
C 8.54338 -2.16689 1.86340
C 8.57094 -3.58890 1.46394
C 8.10584 -3.58300 0.14026
C 7.82847 -2.24915 -0.26395
C 9.03163 -4.66987 2.26730
C 7.79653 -4.40074 -1.01053
O 7.75310 -5.61054 -1.23014
C 7.57792 -3.44768 -2.26364
C 6.31727 -3.87954 -2.96133
O 5.21385 -3.76947 -2.39833
O 6.59617 -4.14595 -4.25036
C 5.37492 -4.56785 -4.98682
H 8.16787 3.43339 -2.95050
H 8.71626 3.93358 3.25116
H 9.17377 -2.30792 3.76977
H 7.84010 -1.24791 -4.48417
H 7.48523 1.54507 -4.58617
H 9.89844 -0.17993 -4.21074
H 10.00922 1.42812 -3.46495
H 9.73000 1.33383 -5.17749
H 5.44114 0.50927 -3.91578
H 5.28283 -1.26485 -3.64626
H 6.36170 -0.86067 -6.25014
H 5.29791 0.47703 -5.99758
H 7.49825 5.40346 -2.61705
H 8.10932 6.56086 -1.48294
H 6.35506 6.26253 -1.50075
H 7.61559 7.62248 0.41563
H 6.78664 8.00124 2.03306
H 5.98272 7.07463 0.73363
H 8.78592 2.28908 5.19806
H 10.39482 -0.19820 4.91242
H 10.87703 3.04520 5.62780
H 11.25621 2.75418 3.89357
H 11.58480 1.54712 5.02965
H 8.86786 0.28533 6.86342
H 8.71673 -1.32484 6.18555
H 7.00057 1.04587 6.05628
H 6.50071 -0.65464 6.37761
H 6.84545 -0.06356 4.69057
H 8.31967 -5.47507 2.14378
H 9.12538 -4.37768 3.31115
H 10.00432 -4.85627 1.80427
H 8.45455 -3.54403 -2.90538
H 5.59370 -5.54092 -5.43735
H 5.13908 -3.73528 -5.67599
H 4.49204 -4.79906 -4.40243

