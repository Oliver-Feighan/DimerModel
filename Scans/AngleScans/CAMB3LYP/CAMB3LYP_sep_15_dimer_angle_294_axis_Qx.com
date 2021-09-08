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
Mg 7.88925 0.70588 0.60660
C 5.41865 -0.08690 3.08585
C 9.27709 2.83160 2.88473
C 9.80523 1.76004 -1.69236
C 5.69635 -0.56577 -1.78500
N 7.32962 1.21916 2.68401
C 6.40802 0.80879 3.54227
C 6.60210 1.32882 4.99191
C 7.46255 2.64408 4.66238
C 8.09490 2.19885 3.34534
C 6.60719 3.87047 4.39100
C 7.40993 0.32269 5.88473
C 7.37146 0.62044 7.38527
H 7.35742 -0.50331 8.31772
N 9.60089 1.70408 0.71833
C 10.01339 2.50183 1.72566
C 11.32991 3.01641 1.31803
C 11.58736 2.58828 -0.08380
C 10.30409 2.00575 -0.42141
C 12.23734 3.76180 2.25311
C 12.78553 2.69796 -0.97743
O 12.70931 2.47153 -2.23670
C 14.02957 3.16873 -0.35115
N 7.68806 0.80820 -1.46013
C 8.69228 1.26139 -2.23090
C 8.42592 0.98582 -3.69155
C 7.14140 0.03363 -3.72184
C 6.81997 -0.00580 -2.21170
C 8.16922 2.21420 -4.56423
C 7.28140 -1.31302 -4.48316
C 8.24995 -2.21832 -3.75671
N 5.99753 -0.33762 0.62951
C 5.24002 -0.75495 -0.48678
C 3.97571 -1.31728 0.03066
C 4.07630 -1.16589 1.42186
C 5.30700 -0.53414 1.74710
C 2.90837 -1.82991 -0.75935
C 3.41656 -1.38720 2.68853
O 2.38879 -1.96270 3.04381
C 4.20056 -0.58381 3.81383
C 4.41085 -1.50707 4.98254
O 5.13780 -2.51012 4.87282
O 3.89876 -0.93979 6.08994
C 4.09281 -1.81292 7.27813
H 9.65594 3.69028 3.43971
H 10.47730 2.00478 -2.51813
H 5.00317 -1.07178 -2.46037
H 5.68663 1.60180 5.50297
H 8.22831 2.81390 5.41480
H 5.55597 3.62705 4.34825
H 6.94444 4.25959 3.43921
H 6.82103 4.64982 5.12804
H 8.48256 0.45636 5.69410
H 7.08940 -0.68561 5.70807
H 6.58716 1.34281 7.64218
H 8.30075 1.12365 7.66314
H 11.71055 4.06945 3.16484
H 12.50020 4.70405 1.76024
H 13.15420 3.21251 2.47742
H 13.86506 4.15123 0.08484
H 14.79353 3.09547 -1.12760
H 14.28140 2.49556 0.47591
H 9.26021 0.46501 -4.14524
H 6.28544 0.51962 -4.21273
H 8.86522 2.29443 -5.42052
H 8.23139 3.09554 -3.93981
H 7.16209 2.28962 -4.97070
H 7.70111 -1.14387 -5.47784
H 6.33921 -1.83431 -4.61557
H 9.22508 -1.87322 -4.11322
H 8.07703 -3.25852 -4.02529
H 8.21443 -2.04538 -2.67483
H 2.58074 -2.74488 -0.28395
H 3.22229 -2.00489 -1.78629
H 2.19539 -1.00463 -0.68302
H 3.59411 0.28260 4.08051
H 3.10461 -1.97701 7.71884
H 4.84550 -1.30008 7.90599
H 4.41837 -2.83174 7.10362

