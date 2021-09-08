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
Mg 7.88072 0.70249 0.60621
C 5.90120 2.67454 2.85823
C 10.01335 0.14213 3.20861
C 9.87466 -0.66105 -1.45107
C 6.26022 2.32476 -2.01828
N 7.85284 1.37487 2.71319
C 7.01839 2.07002 3.47134
C 7.35934 2.05329 4.98562
C 8.92938 1.73448 4.87509
C 8.93697 1.00614 3.53265
C 9.78905 2.98103 4.74704
C 6.59612 0.92040 5.75767
C 6.64093 1.03881 7.28270
H 5.46662 0.62828 8.04750
N 9.40290 -0.53191 0.91896
C 10.14220 -0.64385 2.04243
C 11.15065 -1.67679 1.75918
C 11.05069 -2.04820 0.32155
C 10.09873 -1.05825 -0.14089
C 12.03098 -2.26763 2.82189
C 11.71317 -3.10955 -0.50360
O 11.65795 -3.09932 -1.78410
C 12.50653 -4.10841 0.22739
N 8.19919 0.96216 -1.43180
C 9.09040 0.20407 -2.09436
C 8.95012 0.37704 -3.58795
C 7.61299 1.22881 -3.79810
C 7.24351 1.48833 -2.32110
C 10.11438 1.07677 -4.28885
C 6.53235 0.62265 -4.73483
C 5.94401 -0.62057 -4.10742
N 6.23236 2.08594 0.41375
C 5.73449 2.65922 -0.77686
C 4.68514 3.62432 -0.38973
C 4.66109 3.56191 1.01171
C 5.64117 2.63900 1.46688
C 3.94061 4.44614 -1.28214
C 4.03676 4.07477 2.20999
O 3.08407 4.81747 2.44346
C 4.89858 3.61674 3.46454
C 3.95610 3.06596 4.49937
O 3.31053 2.02655 4.27690
O 4.13453 3.73011 5.65585
C 3.22957 3.21231 6.71639
H 10.86162 0.09132 3.89179
H 10.46167 -1.18438 -2.20937
H 5.64055 2.79844 -2.78268
H 7.20652 2.99660 5.49611
H 9.25448 1.06988 5.67131
H 9.19100 3.87284 4.63230
H 10.40617 2.82418 3.87205
H 10.47842 3.05209 5.59324
H 7.13433 -0.02722 5.62690
H 5.57674 0.85615 5.43001
H 6.98605 2.02729 7.60942
H 7.40047 0.35087 7.66201
H 11.99301 -1.68068 3.74787
H 13.06486 -2.16225 2.47570
H 11.82409 -3.32498 3.00039
H 13.28957 -3.60602 0.79034
H 12.82671 -4.83358 -0.52326
H 11.85855 -4.60135 0.96071
H 8.83703 -0.58153 -4.07933
H 7.82187 2.21098 -4.24702
H 10.56287 0.47197 -5.09950
H 10.85909 1.32805 -3.54534
H 9.87826 2.04936 -4.71723
H 6.98278 0.30988 -5.68004
H 5.73059 1.31356 -4.97320
H 6.66561 -1.39782 -4.37628
H 4.96181 -0.83236 -4.52510
H 5.93505 -0.54310 -3.01400
H 2.91097 4.41131 -0.95179
H 4.04042 4.10761 -2.31131
H 4.43226 5.40859 -1.11709
H 5.43927 4.49163 3.82757
H 2.65792 4.06598 7.09342
H 3.88386 2.68647 7.43690
H 2.43370 2.54225 6.41295

