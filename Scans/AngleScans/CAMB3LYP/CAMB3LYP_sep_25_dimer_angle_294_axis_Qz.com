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
Mg 13.14518 1.19109 0.99733
C 12.55277 2.80962 4.14511
C 12.11689 4.07364 -0.50869
C 13.22539 -0.34676 -1.77683
C 13.02085 -1.84984 2.69528
N 12.44983 3.12328 1.81915
C 12.32491 3.64369 3.03073
C 12.01194 5.16398 3.04859
C 11.34624 5.28957 1.59258
C 12.03566 4.11215 0.90623
C 9.84844 5.03277 1.59147
C 13.30861 6.03985 3.16482
C 13.06192 7.51290 3.49804
H 14.02389 8.20383 4.35227
N 13.20678 1.89205 -0.85831
C 12.76618 3.09247 -1.28975
C 13.01969 3.11707 -2.73865
C 13.49462 1.76946 -3.15470
C 13.30333 1.03061 -1.92279
C 12.88428 4.36137 -3.56739
C 14.02038 1.23130 -4.45091
O 14.15007 -0.02754 -4.65421
C 14.29073 2.21639 -5.50826
N 12.91576 -0.81154 0.48822
C 13.04963 -1.22074 -0.78558
C 13.11680 -2.72696 -0.87045
C 13.24375 -3.23713 0.63991
C 13.13339 -1.88304 1.37453
C 11.92379 -3.40942 -1.53930
C 14.44672 -4.16208 0.97156
C 15.74007 -3.38782 0.85683
N 13.02900 0.56632 3.06233
C 12.96236 -0.75700 3.55067
C 12.77892 -0.67231 5.01409
C 12.72809 0.70534 5.27482
C 12.84959 1.42798 4.05723
C 12.63478 -1.77729 5.89974
C 12.58862 1.71244 6.30194
O 12.53875 1.69844 7.53119
C 12.34200 3.11780 5.60155
C 13.25682 4.12401 6.24432
O 14.49053 4.02448 6.12331
O 12.52388 5.15227 6.70893
C 13.38793 6.18175 7.34564
H 11.59889 4.84499 -1.07929
H 13.35923 -0.92633 -2.69316
H 13.07120 -2.75224 3.30832
H 11.30463 5.47326 3.80874
H 11.62081 6.22503 1.11209
H 9.49756 4.69101 2.55388
H 9.68036 4.27563 0.83680
H 9.31355 5.92335 1.24928
H 13.75985 6.14532 2.16983
H 13.99100 5.61770 3.87669
H 12.03847 7.69001 3.85035
H 13.13203 8.09388 2.57524
H 12.36037 5.15443 -3.01969
H 12.21568 4.11941 -4.40063
H 13.83954 4.70578 -3.96931
H 13.37794 2.76196 -5.73550
H 14.74380 1.65136 -6.32517
H 15.00947 2.95039 -5.12711
H 13.99233 -3.04447 -1.42339
H 12.36977 -3.83174 0.94395
H 12.20929 -4.04071 -2.40187
H 11.22011 -2.64589 -1.84327
H 11.31948 -4.03190 -0.88163
H 14.50811 -4.97795 0.24717
H 14.38768 -4.61067 1.95760
H 15.95745 -3.42789 -0.21472
H 16.52606 -3.86681 1.43720
H 15.59349 -2.33491 1.12430
H 13.24564 -1.56322 6.76664
H 12.93145 -2.70962 5.42379
H 11.55864 -1.73698 6.08834
H 11.29071 3.37297 5.74030
H 13.01158 6.33027 8.36249
H 13.36247 7.05468 6.66659
H 14.42137 5.91889 7.53891

