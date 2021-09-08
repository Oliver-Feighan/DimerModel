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
Mg 7.89124 0.70513 0.59375
C 7.76858 2.41990 -2.55642
C 7.40038 -2.19412 -1.13485
C 7.54448 -0.77901 3.37665
C 7.23025 3.76419 2.12703
N 7.63315 0.29472 -1.56431
C 7.66626 1.01802 -2.67340
C 7.66595 0.18621 -3.98400
C 6.97849 -1.15100 -3.41994
C 7.38801 -1.03631 -1.95304
C 5.46094 -1.12546 -3.49558
C 9.11579 -0.09327 -4.51490
C 9.18412 -0.63045 -5.94616
H 10.30023 -0.21700 -6.79223
N 8.01437 -1.22361 1.04473
C 7.80605 -2.26801 0.21581
C 7.99828 -3.47549 1.03387
C 8.17505 -3.05844 2.45145
C 7.89120 -1.64121 2.34663
C 8.08080 -4.85339 0.44394
C 8.52716 -3.81286 3.69767
O 8.39278 -3.29168 4.86092
C 8.94516 -5.21133 3.52141
N 7.24321 1.34600 2.46235
C 7.23859 0.51105 3.51628
C 6.99833 1.26065 4.80499
C 7.08355 2.81113 4.42246
C 7.27827 2.68471 2.89553
C 5.66558 0.98047 5.49901
C 8.08915 3.68011 5.22639
C 9.50564 3.26201 4.90370
N 7.73879 2.76380 -0.04396
C 7.42836 3.88426 0.75746
C 7.32727 5.04999 -0.14436
C 7.56508 4.51849 -1.42086
C 7.77798 3.11678 -1.32394
C 7.00150 6.37908 0.24745
C 7.67392 4.82271 -2.82951
O 7.67391 5.85583 -3.49764
C 7.68532 3.45671 -3.64211
C 8.80971 3.51883 -4.63930
O 9.99213 3.56414 -4.25662
O 8.31834 3.29596 -5.87189
C 9.39746 3.33660 -6.89447
H 7.02432 -3.12852 -1.55240
H 7.52063 -1.22239 4.37481
H 7.12407 4.76994 2.53918
H 7.06770 0.60208 -4.78582
H 7.41488 -2.03924 -3.86932
H 5.08930 -0.16092 -3.80813
H 5.11271 -1.36002 -2.49820
H 5.10187 -1.93445 -4.13819
H 9.53683 -0.94960 -3.97252
H 9.72879 0.78269 -4.42804
H 8.23279 -0.50829 -6.47814
H 9.33064 -1.71247 -5.90409
H 7.75026 -4.86360 -0.60202
H 7.34285 -5.47033 0.96804
H 9.06802 -5.30462 0.56418
H 8.15178 -5.76563 3.02566
H 9.24084 -5.55533 4.51444
H 9.81199 -5.23560 2.85161
H 7.76827 1.03540 5.53271
H 6.12154 3.32199 4.57534
H 5.78203 0.62151 6.53904
H 5.11741 0.25960 4.90707
H 4.97674 1.82280 5.53336
H 7.95499 3.51996 6.29901
H 7.98269 4.74423 5.04387
H 9.66536 2.40280 5.56198
H 10.20310 4.06642 5.12840
H 9.58482 2.90479 3.87036
H 7.67875 7.03497 -0.28325
H 7.08468 6.50771 1.32462
H 5.96426 6.43429 -0.09366
H 6.71345 3.35569 -4.12672
H 9.10162 4.07933 -7.64181
H 9.51930 2.29341 -7.24182
H 10.36425 3.71846 -6.58797

