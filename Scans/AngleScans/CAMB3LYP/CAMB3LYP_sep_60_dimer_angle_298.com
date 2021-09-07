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
Mg 31.53913 2.84712 2.39509
C 30.95761 4.26125 5.64181
C 30.48882 5.81365 1.07936
C 31.61027 1.48744 -0.47085
C 31.43869 -0.29528 3.89875
N 30.83968 4.72024 3.34084
C 30.71943 5.16270 4.58346
C 30.39956 6.67727 4.69879
C 29.72483 6.89100 3.25729
C 30.41565 5.76258 2.49433
C 28.22825 6.62735 3.24827
C 31.69284 7.55051 4.86282
C 31.44131 8.99841 5.28954
H 32.40505 9.63894 6.18028
N 31.58670 3.66384 0.58698
C 31.13807 4.88685 0.23447
C 31.38302 5.00391 -1.21137
C 31.86173 3.68755 -1.71412
C 31.68102 2.87164 -0.53017
C 31.23706 6.29724 -1.95932
C 32.38241 3.23470 -3.04454
O 32.51670 1.99181 -3.32746
C 32.64205 4.28575 -4.03921
N 31.31598 0.87941 1.76209
C 31.44431 0.55191 0.46433
C 31.51792 -0.94562 0.28436
C 31.65601 -1.54926 1.75886
C 31.54369 -0.24469 2.57792
C 30.32420 -1.59052 -0.41956
C 32.86514 -2.48728 2.02492
C 34.15422 -1.70092 1.95207
N 31.43785 2.09297 4.41724
C 31.38015 0.74122 4.82159
C 31.20485 0.73266 6.28843
C 31.14919 2.09089 6.63571
C 31.26027 2.88937 5.46542
C 31.07096 -0.42659 7.10348
C 31.01106 3.03061 7.72500
O 30.96842 2.93898 8.95118
C 30.75390 4.47604 7.11591
C 31.66780 5.44431 7.81574
O 32.90124 5.35872 7.68189
O 30.93285 6.43762 8.34824
C 31.79585 7.42924 9.04377
H 29.96396 6.61682 0.56136
H 31.74144 0.96739 -1.42260
H 31.49676 -1.23423 4.45342
H 29.69527 6.93455 5.48082
H 29.99228 7.85622 2.83518
H 27.88455 6.22392 4.18915
H 28.05926 5.91841 2.44834
H 27.68727 7.53504 2.96583
H 32.13778 7.72066 3.87397
H 32.38130 7.08776 5.54290
H 30.41913 9.14790 5.65796
H 31.50337 9.63669 4.40481
H 30.71270 7.05162 -1.35986
H 30.56474 6.10491 -2.80244
H 32.18836 6.67101 -2.34401
H 31.72545 4.84000 -4.22658
H 33.09296 3.77554 -4.89258
H 33.35961 4.99785 -3.61654
H 32.39167 -1.22334 -0.29231
H 30.78657 -2.16617 2.02965
H 30.60758 -2.16481 -1.32175
H 29.61525 -0.81287 -0.67094
H 29.72659 -2.25617 0.20091
H 32.92607 -3.25560 1.25024
H 32.81391 -2.99736 2.98106
H 34.36554 -1.67235 0.87893
H 34.94577 -2.21160 2.49676
H 34.00435 -0.66768 2.28614
H 31.68586 -0.26451 7.97877
H 31.36914 -1.32562 6.56812
H 29.99576 -0.40359 7.30020
H 29.70227 4.71675 7.27627
H 31.42474 7.51156 10.07002
H 31.76241 8.34305 8.42122
H 32.83159 7.15985 9.21439

