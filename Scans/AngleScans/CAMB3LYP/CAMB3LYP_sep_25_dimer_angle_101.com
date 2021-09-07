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
Mg 13.14619 1.17633 1.00151
C 13.35470 -1.13742 -1.73384
C 12.35085 -1.47789 2.99077
C 12.51480 3.23753 3.32965
C 12.82711 3.72176 -1.35745
N 12.91427 -0.98026 0.56841
C 13.07856 -1.75260 -0.49500
C 13.01844 -3.27530 -0.19949
C 12.14650 -3.22448 1.14811
C 12.51522 -1.82118 1.62503
C 10.64824 -3.26545 0.89739
C 14.43863 -3.89070 0.05813
C 14.48840 -5.42013 0.04711
H 15.67733 -6.07397 -0.49271
N 13.00623 0.89678 2.96131
C 12.69580 -0.25115 3.59941
C 12.70337 0.07163 5.03466
C 12.88171 1.54153 5.18394
C 12.78487 1.96071 3.80010
C 12.63023 -0.97817 6.10517
C 13.09042 2.42527 6.37627
O 12.98262 3.69998 6.29672
C 13.33256 1.74911 7.65910
N 12.51960 3.15868 0.99810
C 12.37125 3.83881 2.14847
C 12.18320 5.31456 1.88929
C 12.47900 5.51819 0.33089
C 12.70984 4.04662 -0.07726
C 10.80191 5.87782 2.22270
C 13.55820 6.56867 -0.04923
C 14.92028 6.09759 0.40706
N 13.28026 1.32927 -1.14946
C 13.08719 2.48604 -1.93603
C 13.16758 2.06886 -3.35092
C 13.38060 0.68347 -3.28905
C 13.40937 0.26351 -1.93174
C 13.00040 2.91958 -4.47984
C 13.57707 -0.52081 -4.06351
O 13.73248 -0.76953 -5.25839
C 13.44281 -1.77289 -3.09363
C 14.60047 -2.69489 -3.36244
O 15.76474 -2.33760 -3.11057
O 14.12855 -3.91786 -3.66577
C 15.23949 -4.87069 -3.92985
H 11.87364 -2.20033 3.65352
H 12.39932 4.00769 4.09580
H 12.83562 4.47144 -2.15154
H 12.50685 -3.86411 -0.95138
H 12.48161 -3.97074 1.86377
H 10.41397 -3.20261 -0.15490
H 10.23791 -2.41608 1.42772
H 10.21183 -4.15248 1.36542
H 14.72755 -3.70082 1.09986
H 15.15492 -3.50082 -0.63868
H 13.57982 -5.85903 -0.38279
H 14.49410 -5.77460 1.08072
H 12.33825 -1.95236 5.69395
H 11.80154 -0.70315 6.76662
H 13.54694 -1.04265 6.69522
H 12.49297 1.09698 7.88772
H 13.54666 2.54620 8.37373
H 14.21234 1.10444 7.55420
H 12.89191 5.90101 2.46106
H 11.58522 5.85849 -0.21243
H 10.83479 6.71569 2.94448
H 10.18751 5.07244 2.60273
H 10.22538 6.22399 1.36656
H 13.36689 7.51224 0.46771
H 13.59489 6.78520 -1.11166
H 14.94572 6.39796 1.45883
H 15.70630 6.58893 -0.16286
H 14.98971 5.00407 0.37477
H 13.77412 2.65332 -5.18762
H 13.06119 3.96904 -4.19929
H 11.99139 2.63549 -4.79041
H 12.48373 -2.24816 -3.30308
H 15.06742 -5.29442 -4.92410
H 15.23944 -5.57375 -3.07576
H 16.23566 -4.45970 -4.04462

