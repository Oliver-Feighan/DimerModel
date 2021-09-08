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
Mg 13.14649 1.17580 0.99816
C 13.28480 0.25919 -2.46880
C 12.50631 -2.10327 1.68551
C 12.57327 2.02037 4.00241
C 12.65944 4.46837 -0.03500
N 12.94643 -0.59806 -0.30785
C 13.08028 -0.83607 -1.60392
C 13.07000 -2.34055 -1.98600
C 12.25960 -2.90474 -0.71959
C 12.61621 -1.82469 0.29972
C 10.75284 -2.89640 -0.91694
C 14.51472 -2.94819 -2.06058
C 14.60033 -4.32360 -2.72591
H 15.77857 -4.63561 -3.53019
N 13.10294 0.08217 2.65365
C 12.84941 -1.23990 2.74910
C 12.91488 -1.55958 4.18356
C 13.06481 -0.28753 4.94126
C 12.89488 0.67708 3.87320
C 12.91574 -2.96753 4.70424
C 13.30670 0.01152 6.38980
O 13.16506 1.19295 6.86622
C 13.62324 -1.13632 7.25221
N 12.47339 2.94299 1.86205
C 12.36162 3.06126 3.19676
C 12.12684 4.49764 3.59934
C 12.34618 5.35797 2.26919
C 12.59309 4.21153 1.26415
C 10.74921 4.80787 4.18439
C 13.38157 6.51359 2.34122
C 14.77391 5.94932 2.50965
N 13.17844 2.23624 -0.88400
C 12.92214 3.60891 -1.09411
C 12.94770 3.83821 -2.55329
C 13.19624 2.56861 -3.09627
C 13.29700 1.61180 -2.05040
C 12.70890 5.08128 -3.20434
C 13.38572 1.81843 -4.31693
O 13.49224 2.10928 -5.50765
C 13.32577 0.26806 -3.97175
C 14.49156 -0.40308 -4.64483
O 15.65728 -0.13975 -4.30080
O 14.03551 -1.39837 -5.42699
C 15.15559 -2.10130 -6.10753
H 12.07722 -3.05819 1.99032
H 12.47463 2.38510 5.02746
H 12.61379 5.48464 -0.43222
H 12.53874 -2.57320 -2.90110
H 12.64473 -3.87039 -0.40247
H 10.46930 -2.40076 -1.83350
H 10.34707 -2.37169 -0.06191
H 10.35951 -3.91540 -0.85970
H 14.84634 -3.20870 -1.04722
H 15.18893 -2.26948 -2.54578
H 13.68376 -4.57430 -3.27364
H 14.66170 -5.08421 -1.94378
H 12.62856 -3.68464 3.92534
H 12.11184 -3.03491 5.44523
H 13.85974 -3.23965 5.18121
H 12.81073 -1.85755 7.20613
H 13.85074 -0.71177 8.23194
H 14.51240 -1.63809 6.85443
H 12.84678 4.81302 4.34455
H 11.42066 5.86038 1.95164
H 10.79510 5.25879 5.19367
H 10.17217 3.89286 4.20259
H 10.12609 5.46192 3.57683
H 13.19168 7.13813 3.21760
H 13.36450 7.16358 1.47276
H 14.84023 5.77356 3.58752
H 15.52114 6.66849 2.18030
H 14.86780 4.97771 2.01084
H 13.45559 5.17407 -3.98174
H 12.75745 5.91261 -2.50413
H 11.69380 4.91550 -3.57496
H 12.36971 -0.11154 -4.33427
H 14.94840 -2.06756 -7.18160
H 15.21131 -3.10066 -5.63644
H 16.13541 -1.64009 -6.06656

