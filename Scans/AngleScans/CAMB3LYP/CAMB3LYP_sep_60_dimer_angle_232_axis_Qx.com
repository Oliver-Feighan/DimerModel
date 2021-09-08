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
Mg 31.53115 2.83142 2.40281
C 29.34570 4.56504 4.66057
C 33.74716 2.57283 4.98293
C 33.65671 1.70551 0.33357
C 29.68953 4.20748 -0.21648
N 31.44334 3.51898 4.50324
C 30.53697 4.11253 5.26518
C 30.89478 4.15599 6.77503
C 32.49104 4.03548 6.64744
C 32.57458 3.29858 5.31226
C 33.18578 5.37851 6.49603
C 30.28887 2.94516 7.56798
C 30.33605 3.08564 9.09107
H 29.23150 2.53974 9.87494
N 33.19984 1.80145 2.70809
C 33.96026 1.79599 3.82300
C 35.08713 0.89460 3.53699
C 35.01805 0.49720 2.10450
C 33.94400 1.35453 1.64461
C 36.04691 0.43106 4.59405
C 35.79904 -0.48196 1.28140
O 35.72820 -0.49335 0.00168
C 36.71997 -1.36497 2.01194
N 31.79095 3.10575 0.35840
C 32.76262 2.45798 -0.30807
C 32.58450 2.59493 -1.80136
C 31.14863 3.26972 -2.00266
C 30.76654 3.49766 -0.52366
C 33.64347 3.42719 -2.52407
C 30.14195 2.52210 -2.91937
C 29.72173 1.22211 -2.27221
N 29.71995 3.99481 2.21795
C 29.14031 4.48747 1.02825
C 27.98257 5.31761 1.41932
C 27.98272 5.26866 2.82150
C 29.07616 4.48129 3.27310
C 27.13045 6.02927 0.52851
C 27.31278 5.71275 4.02261
O 26.27709 6.33265 4.26098
C 28.23972 5.38083 5.27041
C 27.38590 4.72797 6.32268
O 26.87349 3.61331 6.11881
O 27.49282 5.42241 7.47015
C 26.67237 4.80728 8.54731
H 34.60292 2.63667 5.65560
H 34.29602 1.25137 -0.42702
H 29.00649 4.59090 -0.97747
H 30.63057 5.07842 7.27811
H 32.90622 3.42605 7.44593
H 32.47916 6.18685 6.38022
H 33.80758 5.29037 5.61476
H 33.87048 5.54516 7.33254
H 30.94031 2.07114 7.43965
H 29.28192 2.74977 7.25413
H 30.55800 4.11328 9.40355
H 31.18031 2.50283 9.46733
H 35.94620 1.01913 5.51461
H 37.05531 0.66139 4.23355
H 35.97656 -0.64179 4.78563
H 37.44014 -0.76190 2.55977
H 37.12002 -2.05274 1.26441
H 36.14755 -1.92693 2.75836
H 32.58706 1.62421 -2.28176
H 31.22728 4.26516 -2.46392
H 34.15500 2.87424 -3.33442
H 34.35924 3.77840 -1.79272
H 33.28212 4.35751 -2.95895
H 30.61718 2.25757 -3.86718
H 29.25706 3.10418 -3.15421
H 30.53209 0.53855 -2.54265
H 28.76915 0.88400 -2.67509
H 29.71573 1.31032 -1.17958
H 26.11723 5.86927 0.87242
H 27.26012 5.69422 -0.49847
H 27.49918 7.04763 0.67771
H 28.67039 6.32074 5.61780
H 26.00240 5.58671 8.92324
H 27.39580 4.37595 9.26450
H 25.96353 4.03923 8.26078

