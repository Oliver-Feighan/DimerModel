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
Mg 8.94284 0.80909 0.68261
C 9.69574 -1.33788 3.45799
C 8.60749 3.35093 2.93232
C 7.83919 2.61974 -1.67767
C 8.22026 -2.05880 -1.16098
N 9.15041 0.87389 2.88374
C 9.49997 0.00991 3.82485
C 9.72643 0.64065 5.22501
C 8.81322 1.94990 5.04951
C 8.88462 2.09593 3.53087
C 7.36154 1.72789 5.44027
C 11.22893 1.01760 5.47485
C 11.57547 1.35023 6.92778
H 12.88751 0.96083 7.43707
N 8.78930 2.78272 0.54202
C 8.68516 3.66522 1.55772
C 8.57806 4.99105 0.92934
C 8.46171 4.80454 -0.54253
C 8.33520 3.36279 -0.61635
C 8.67185 6.27075 1.70845
C 8.45148 5.76583 -1.69231
O 8.10084 5.40272 -2.87046
C 8.77402 7.16566 -1.37938
N 7.94333 0.36564 -1.08563
C 7.62418 1.33493 -1.96115
C 7.16228 0.75231 -3.27553
C 7.46886 -0.81450 -3.18235
C 7.99597 -0.88342 -1.73225
C 5.68672 0.96262 -3.61449
C 8.33623 -1.43043 -4.31416
C 9.74654 -0.89183 -4.23476
N 9.12211 -1.32230 0.98771
C 8.73633 -2.34705 0.09600
C 8.94728 -3.63260 0.79262
C 9.42321 -3.26276 2.05958
C 9.48410 -1.84608 2.15350
C 8.65872 -4.92235 0.26403
C 9.87787 -3.74838 3.34261
O 10.12167 -4.85808 3.81470
C 9.95471 -2.52074 4.34928
C 11.27846 -2.58551 5.06061
O 12.34119 -2.43057 4.43354
O 11.06457 -2.60204 6.38890
C 12.34835 -2.65486 7.13785
H 8.25625 4.16326 3.56923
H 7.54862 3.19415 -2.56034
H 8.11149 -3.00112 -1.70214
H 9.36654 0.04485 6.05523
H 9.26098 2.81203 5.53694
H 7.15764 0.69049 5.66020
H 6.77505 2.05729 4.59258
H 7.08933 2.38713 6.26951
H 11.43765 1.98743 5.00515
H 11.88041 0.24476 5.11581
H 10.78567 1.03793 7.62179
H 11.61269 2.43697 7.03560
H 8.58988 6.09160 2.78761
H 7.78220 6.86066 1.46269
H 9.56171 6.85216 1.45791
H 8.06955 7.54234 -0.64162
H 8.80314 7.68093 -2.34142
H 9.76555 7.20068 -0.91434
H 7.72235 1.17142 -4.10243
H 6.54617 -1.41259 -3.21008
H 5.53002 1.47727 -4.58127
H 5.22733 1.51931 -2.80852
H 5.08516 0.05563 -3.64103
H 7.94660 -1.13717 -5.29208
H 8.36841 -2.51462 -4.29225
H 9.67514 0.06524 -4.76012
H 10.44207 -1.56445 -4.73262
H 10.02837 -0.67809 -3.19719
H 9.49446 -5.55889 0.52257
H 8.50416 -4.88542 -0.81231
H 7.73607 -5.15230 0.80357
H 9.11453 -2.60961 5.03904
H 12.29793 -3.52688 7.79720
H 12.45419 -1.66463 7.61951
H 13.24909 -2.86777 6.57412

