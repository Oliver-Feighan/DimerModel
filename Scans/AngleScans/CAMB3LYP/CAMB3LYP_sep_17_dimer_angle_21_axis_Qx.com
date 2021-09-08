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
Mg 8.93850 0.81148 0.68522
C 9.63788 -1.22782 3.55423
C 6.93023 2.69433 2.69926
C 7.95111 2.35192 -1.90699
C 10.04845 -1.83437 -1.29295
N 8.43089 0.65394 2.83226
C 8.80510 -0.12489 3.83613
C 8.32103 0.34894 5.23281
C 7.04905 1.21128 4.76666
C 7.50556 1.57721 3.35586
C 5.77359 0.39222 4.65858
C 9.38455 1.24247 5.96258
C 9.11064 1.48407 7.44864
H 10.24762 1.57005 8.36097
N 8.03677 2.56982 0.50164
C 7.27979 3.19661 1.42652
C 6.83063 4.44837 0.79763
C 7.24827 4.43193 -0.63066
C 7.75444 3.07858 -0.74173
C 6.15158 5.54605 1.56406
C 7.18869 5.45945 -1.72008
O 7.38925 5.15136 -2.94796
C 6.79685 6.81846 -1.31876
N 8.79038 0.25778 -1.31323
C 8.37938 1.13556 -2.24517
C 8.61527 0.60255 -3.63824
C 9.49565 -0.71911 -3.44920
C 9.54317 -0.77299 -1.90640
C 7.36014 0.26682 -4.44343
C 10.83481 -0.79075 -4.23290
C 11.80356 0.23733 -3.69443
N 9.87797 -1.10287 1.03265
C 10.23455 -2.06638 0.06386
C 10.73609 -3.25031 0.79133
C 10.61400 -2.89940 2.14427
C 10.05563 -1.59711 2.25275
C 11.17654 -4.46664 0.19748
C 10.82061 -3.33869 3.50553
O 11.34539 -4.31924 4.03172
C 10.08018 -2.32577 4.48131
C 11.03318 -1.96663 5.58825
O 12.07222 -1.32789 5.34533
O 10.45778 -2.24416 6.77243
C 11.36019 -1.89475 7.90176
H 6.10292 3.21232 3.18519
H 7.72926 2.88182 -2.83622
H 10.50455 -2.66249 -1.83964
H 8.00746 -0.44607 5.89870
H 6.93385 2.10525 5.37405
H 5.95880 -0.66229 4.79944
H 5.38728 0.57911 3.66515
H 5.01962 0.77539 5.35215
H 9.30612 2.27099 5.58743
H 10.37081 0.83957 5.83777
H 8.35120 0.79888 7.84452
H 8.66342 2.47455 7.56251
H 5.83128 5.20518 2.55634
H 5.21677 5.77135 1.03943
H 6.75357 6.45554 1.61985
H 5.81440 6.78851 -0.85361
H 6.89832 7.42936 -2.21785
H 7.49663 7.17275 -0.55365
H 9.17285 1.31381 -4.23517
H 8.95564 -1.61896 -3.77858
H 7.30332 0.80561 -5.40808
H 6.49419 0.48695 -3.83334
H 7.22326 -0.79001 -4.66608
H 10.67186 -0.54273 -5.28461
H 11.30185 -1.76946 -4.20145
H 11.50593 1.15292 -4.21423
H 12.82813 -0.03835 -3.93589
H 11.64732 0.40152 -2.62192
H 12.08029 -4.75712 0.71647
H 11.35086 -4.34768 -0.86987
H 10.31516 -5.10753 0.40380
H 9.18676 -2.82631 4.85659
H 11.47782 -2.79823 8.50804
H 10.89970 -1.01607 8.39127
H 12.39209 -1.66278 7.66547

