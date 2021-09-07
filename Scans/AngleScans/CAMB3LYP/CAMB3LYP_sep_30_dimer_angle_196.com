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
Mg 15.77424 1.41777 1.19086
C 15.46697 4.30495 -0.91836
C 15.27406 -0.48572 -1.59502
C 15.59046 -1.12372 3.08149
C 15.11138 3.50165 3.90449
N 15.42764 1.95603 -0.92631
C 15.39280 3.08372 -1.62017
C 15.35240 2.89056 -3.15985
C 14.72347 1.41351 -3.20048
C 15.19675 0.90879 -1.83889
C 13.20386 1.40653 -3.21062
C 16.78271 2.92376 -3.80416
C 16.79832 3.05111 -5.32905
H 17.86444 3.83114 -5.95148
N 15.96373 -0.51239 0.76995
C 15.74268 -1.11149 -0.41895
C 16.00077 -2.54345 -0.20204
C 16.23217 -2.76346 1.25150
C 15.91017 -1.44969 1.77150
C 16.08903 -3.53385 -1.32665
C 16.65872 -3.96191 2.04395
O 16.56528 -3.99220 3.32191
C 17.10140 -5.13342 1.27383
N 15.19722 1.17395 3.17327
C 15.26064 -0.02996 3.76868
C 15.06176 0.08836 5.26094
C 15.09247 1.65606 5.57549
C 15.22016 2.20057 4.13584
C 13.76916 -0.51532 5.80971
C 16.11279 2.14001 6.64200
C 17.52260 1.95783 6.12758
N 15.54384 3.54372 1.49953
C 15.24378 4.20199 2.71227
C 15.07386 5.63566 2.39899
C 15.26568 5.70913 1.01108
C 15.51608 4.40983 0.49278
C 14.73479 6.65640 3.33132
C 15.30275 6.58887 -0.13495
O 15.24761 7.80726 -0.29693
C 15.30951 5.70120 -1.45343
C 16.38536 6.22848 -2.36285
O 17.58262 6.15487 -2.03457
O 15.84362 6.53222 -3.55643
C 16.87360 7.04905 -4.49656
H 14.90165 -1.16762 -2.36001
H 15.62281 -1.95079 3.79451
H 15.00024 4.23053 4.71012
H 14.70840 3.58355 -3.68786
H 15.15988 0.82023 -3.99980
H 12.79548 2.39614 -3.06908
H 12.90726 0.75520 -2.39893
H 12.83516 0.93459 -4.12586
H 17.24836 1.93602 -3.69333
H 17.37801 3.70352 -3.37031
H 15.82103 3.34918 -5.72781
H 16.97233 2.06136 -5.75816
H 15.71137 -3.11087 -2.26572
H 15.39070 -4.34503 -1.09388
H 17.09117 -3.95236 -1.44167
H 16.29962 -5.45564 0.61369
H 17.45026 -5.85539 2.01474
H 17.93705 -4.83433 0.63130
H 15.86929 -0.39374 5.79813
H 14.12656 2.01309 5.96199
H 13.94154 -1.27818 6.59227
H 13.21183 -0.93700 4.98381
H 13.06275 0.20298 6.22230
H 16.03164 1.53263 7.54680
H 15.97279 3.17507 6.93531
H 17.73265 0.90737 6.35004
H 18.21026 2.61765 6.65275
H 17.56296 2.07867 5.03881
H 15.37128 7.50310 3.11113
H 14.86403 6.31698 4.35692
H 13.68203 6.80910 3.07905
H 14.31920 5.77656 -1.90417
H 16.52632 8.02662 -4.84503
H 17.00424 6.25936 -5.26015
H 17.84403 7.30319 -4.08652

