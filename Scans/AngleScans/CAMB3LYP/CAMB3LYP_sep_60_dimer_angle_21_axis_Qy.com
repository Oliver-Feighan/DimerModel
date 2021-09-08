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
Mg 31.54117 2.84356 2.39926
C 32.29408 0.69660 5.17464
C 31.20582 5.38541 4.64897
C 30.43752 4.65421 0.03898
C 30.81859 -0.02432 0.55568
N 31.74874 2.90836 4.60039
C 32.09830 2.04438 5.54150
C 32.32476 2.67513 6.94166
C 31.41155 3.98437 6.76616
C 31.48295 4.13041 5.24752
C 29.95988 3.76237 7.15692
C 33.82726 3.05208 7.19150
C 34.17380 3.38470 8.64444
H 35.48585 2.99531 9.15372
N 31.38763 4.81720 2.25867
C 31.28350 5.69969 3.27437
C 31.17640 7.02552 2.64599
C 31.06004 6.83901 1.17413
C 30.93353 5.39727 1.10030
C 31.27018 8.30522 3.42510
C 31.04981 7.80030 0.02434
O 30.69917 7.43720 -1.15381
C 31.37235 9.20014 0.33727
N 30.54166 2.40011 0.63102
C 30.22252 3.36941 -0.24450
C 29.76061 2.78679 -1.55888
C 30.06719 1.21997 -1.46570
C 30.59430 1.15105 -0.01559
C 28.28505 2.99709 -1.89784
C 30.93456 0.60404 -2.59750
C 32.34488 1.14264 -2.51811
N 31.72044 0.71217 2.70436
C 31.33466 -0.31257 1.81265
C 31.54561 -1.59812 2.50927
C 32.02154 -1.22828 3.77623
C 32.08243 0.18839 3.87016
C 31.25705 -2.88787 1.98068
C 32.47620 -1.71391 5.05926
O 32.72000 -2.82360 5.53135
C 32.55304 -0.48627 6.06593
C 33.87680 -0.55103 6.77726
O 34.93952 -0.39610 6.15019
O 33.66290 -0.56756 8.10555
C 34.94669 -0.62039 8.85450
H 30.85459 6.19773 5.28588
H 30.14695 5.22862 -0.84369
H 30.70982 -0.96665 0.01451
H 31.96487 2.07932 7.77188
H 31.85931 4.84651 7.25359
H 29.75598 2.72497 7.37686
H 29.37338 4.09176 6.30924
H 29.68766 4.42160 7.98616
H 34.03598 4.02190 6.72180
H 34.47874 2.27923 6.83246
H 33.38400 3.07241 9.33844
H 34.21102 4.47145 8.75225
H 31.18821 8.12607 4.50426
H 30.38053 8.89513 3.17934
H 32.16004 8.88663 3.17456
H 30.66788 9.57681 1.07503
H 31.40147 9.71541 -0.62477
H 32.36388 9.23515 0.80231
H 30.32068 3.20589 -2.38578
H 29.14450 0.62189 -1.49343
H 28.12836 3.51174 -2.86462
H 27.82566 3.55378 -1.09187
H 27.68349 2.09011 -1.92438
H 30.54493 0.89731 -3.57543
H 30.96674 -0.48015 -2.57560
H 32.27347 2.09971 -3.04347
H 33.04040 0.47002 -3.01597
H 32.62670 1.35638 -1.48054
H 32.09280 -3.52442 2.23922
H 31.10249 -2.85095 0.90434
H 30.33440 -3.11782 2.52022
H 31.71286 -0.57514 6.75569
H 34.89626 -1.49240 9.51385
H 35.05253 0.36984 9.33616
H 35.84742 -0.83330 8.29077

