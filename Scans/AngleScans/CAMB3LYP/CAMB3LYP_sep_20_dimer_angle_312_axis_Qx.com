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
Mg 10.51805 0.94477 0.80608
C 8.47988 -0.49554 3.38492
C 11.23744 3.46487 2.98923
C 11.95277 2.46769 -1.57924
C 8.76590 -1.01586 -1.48119
N 9.88097 1.32407 2.88995
C 9.15448 0.67577 3.78787
C 9.21538 1.27499 5.21853
C 9.61707 2.78092 4.83160
C 10.32188 2.51233 3.50353
C 8.41704 3.67387 4.56344
C 10.31809 0.59598 6.10453
C 10.22817 0.91361 7.59878
H 10.58723 -0.13000 8.55504
N 11.83841 2.42606 0.83768
C 12.00937 3.34306 1.81296
C 13.09060 4.22644 1.34954
C 13.43160 3.85562 -0.05079
C 12.38371 2.89483 -0.33168
C 13.74631 5.24440 2.23679
C 14.51321 4.30248 -0.98708
O 14.47828 4.02468 -2.23785
C 15.56594 5.15387 -0.41401
N 10.24153 0.91579 -1.25427
C 11.03555 1.63305 -2.06852
C 10.82992 1.24352 -3.51295
C 9.90335 -0.05952 -3.47880
C 9.64923 -0.14954 -1.95804
C 10.18266 2.30475 -4.40246
C 10.43396 -1.31988 -4.21541
C 11.65386 -1.85854 -3.50333
N 9.04406 -0.63108 0.91613
C 8.42447 -1.29648 -0.16435
C 7.41052 -2.20571 0.40803
C 7.49529 -1.98760 1.79140
C 8.47764 -0.99666 2.06068
C 6.53446 -3.04735 -0.33378
C 6.96968 -2.36262 3.08438
O 6.18040 -3.21629 3.48691
C 7.49504 -1.32178 4.16468
C 8.01130 -2.09816 5.34504
O 9.01019 -2.83031 5.23201
O 7.37752 -1.68288 6.45690
C 7.86332 -2.41602 7.65620
H 11.34582 4.41538 3.51237
H 12.49427 2.88245 -2.43260
H 8.24637 -1.73201 -2.12142
H 8.27397 1.26736 5.75468
H 10.31175 3.20232 5.55354
H 7.49226 3.11627 4.56195
H 8.59232 4.11847 3.59254
H 8.39789 4.50360 5.27585
H 11.29118 1.04863 5.87445
H 10.32131 -0.46699 5.96061
H 9.26557 1.36581 7.86686
H 10.96266 1.68775 7.83380
H 13.17395 5.40228 3.15922
H 13.69138 6.20601 1.71511
H 14.79374 5.01261 2.44129
H 15.11643 6.05044 0.00619
H 16.29459 5.29629 -1.21442
H 16.03533 4.61747 0.41818
H 11.77246 0.99226 -3.98367
H 8.92650 0.12267 -3.95039
H 10.79709 2.56958 -5.28361
H 9.98483 3.18106 -3.79958
H 9.19156 2.05257 -4.77572
H 10.75465 -1.06021 -5.22728
H 9.69658 -2.11069 -4.30427
H 12.46447 -1.24018 -3.90034
H 11.80488 -2.90905 -3.74345
H 11.59456 -1.67161 -2.42481
H 6.51893 -4.00358 0.17206
H 6.86041 -3.14847 -1.36689
H 5.60307 -2.48088 -0.25079
H 6.65709 -0.67739 4.43334
H 6.98639 -2.86376 8.13391
H 8.43606 -1.67646 8.24677
H 8.48393 -3.28924 7.49252

