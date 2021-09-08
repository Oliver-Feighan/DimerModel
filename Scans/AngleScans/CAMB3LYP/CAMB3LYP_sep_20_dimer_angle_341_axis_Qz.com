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
Mg 10.51715 0.95297 0.80372
C 10.16642 -0.28405 4.15422
C 9.28074 3.92533 1.93109
C 10.38931 1.93017 -2.21226
C 10.62750 -2.36364 -0.26102
N 9.84671 1.61278 2.80557
C 9.81712 1.07546 4.01584
C 9.46627 2.07163 5.15343
C 8.66943 3.16034 4.28243
C 9.32584 2.91554 2.92528
C 7.18502 2.86048 4.15607
C 10.74413 2.69250 5.81947
C 10.48886 3.43267 7.13428
H 11.50471 3.36267 8.18100
N 10.39307 2.78531 0.05159
C 9.88342 3.87971 0.65476
C 10.00509 4.97165 -0.32352
C 10.47680 4.39633 -1.61240
C 10.41679 2.98099 -1.30696
C 9.76245 6.41056 0.02876
C 10.89835 5.01678 -2.90991
O 11.04304 4.31817 -3.97473
C 11.04612 6.47951 -2.91971
N 10.29685 -0.05992 -0.99847
C 10.32679 0.59915 -2.17007
C 10.42652 -0.35933 -3.33272
C 10.70220 -1.79474 -2.68375
C 10.62188 -1.41611 -1.18856
C 9.19706 -0.43781 -4.23752
C 11.95450 -2.56493 -3.18496
C 13.21105 -1.84555 -2.75026
N 10.60364 -0.98434 1.75630
C 10.61669 -2.24643 1.12300
C 10.56315 -3.26862 2.18822
C 10.49893 -2.52424 3.37577
C 10.49117 -1.13621 3.07112
C 10.52884 -4.67750 1.98823
C 10.42521 -2.59714 4.81728
O 10.48628 -3.50499 5.64542
C 10.07903 -1.15124 5.37935
C 11.02049 -0.85665 6.51487
O 12.24059 -0.73283 6.30840
O 10.30485 -0.55581 7.61384
C 11.19459 -0.24524 8.76437
H 8.69306 4.82222 2.12884
H 10.45588 2.21320 -3.26533
H 10.75691 -3.41992 -0.50635
H 8.82200 1.66996 5.92632
H 8.87452 4.17002 4.62861
H 6.93126 1.89822 4.57536
H 6.97027 2.87963 3.09562
H 6.59783 3.67103 4.59716
H 11.10123 3.52595 5.20094
H 11.49876 1.94448 5.96650
H 9.49678 3.21110 7.54614
H 10.46021 4.50548 6.92867
H 9.26869 6.50832 1.00346
H 9.02850 6.79571 -0.68737
H 10.66816 7.01769 -0.03241
H 10.10229 6.93972 -2.63714
H 11.43905 6.72692 -3.90777
H 11.77619 6.76252 -2.15321
H 11.25709 -0.09987 -3.97766
H 9.87534 -2.49426 -2.87563
H 9.42084 -0.21637 -5.29816
H 8.44867 0.24424 -3.85637
H 8.67116 -1.39088 -4.22042
H 11.97252 -2.58854 -4.27733
H 11.99629 -3.59264 -2.83989
H 13.33237 -1.07516 -3.51764
H 14.05856 -2.52793 -2.73858
H 13.06083 -1.33529 -1.79176
H 11.20908 -5.11097 2.70925
H 10.80653 -4.94139 0.96989
H 9.47328 -4.87759 2.19080
H 9.03797 -1.16665 5.70420
H 10.90718 -0.91532 9.58051
H 11.08483 0.84108 8.94218
H 12.24790 -0.47808 8.66054

