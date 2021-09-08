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
Mg 8.92730 0.80817 0.68572
C 6.55906 -1.88355 0.84439
C 6.36423 2.80344 -0.35550
C 11.06413 3.08041 0.10407
C 11.44789 -1.59773 0.62240
N 6.78330 0.39759 0.32945
C 5.98409 -0.64749 0.48244
C 4.47370 -0.33260 0.31175
C 4.60373 0.97965 -0.60475
C 5.98800 1.45011 -0.16303
C 4.65584 0.67220 -2.09206
C 3.77188 -0.01180 1.67797
C 2.24248 0.00114 1.62787
H 1.50386 -0.48316 2.79078
N 8.69054 2.76497 0.45460
C 7.57610 3.40715 0.04604
C 7.92189 4.83683 0.02108
C 9.37943 4.97166 0.28957
C 9.78059 3.57923 0.27218
C 6.89727 5.91876 -0.16062
C 10.26777 6.15827 0.51099
O 11.54526 6.05628 0.48992
C 9.59907 7.45906 0.66021
N 10.94579 0.75215 0.19084
C 11.65343 1.88605 0.04445
C 13.13369 1.59861 -0.03675
C 13.29136 0.04995 0.32980
C 11.80123 -0.32780 0.47888
C 13.79093 1.87058 -1.38969
C 14.26304 -0.30354 1.48891
C 13.71234 0.21146 2.79922
N 9.03489 -1.33805 0.90983
C 10.18831 -2.14803 0.82245
C 9.74297 -3.55264 0.92819
C 8.34792 -3.46268 1.04758
C 7.94995 -2.09935 0.99803
C 10.58357 -4.69928 0.85939
C 7.12053 -4.21160 1.19358
O 6.84212 -5.39592 1.37693
C 5.89635 -3.22966 0.94147
C 4.89677 -3.44050 2.04551
O 5.18197 -3.14928 3.22044
O 3.70206 -3.74408 1.50619
C 2.67483 -3.95155 2.56145
H 5.68556 3.45772 -0.90343
H 11.85298 3.83039 0.01071
H 12.18193 -2.40162 0.70965
H 3.90672 -1.09503 -0.20883
H 3.84948 1.71832 -0.34621
H 4.71597 -0.38911 -2.28221
H 5.53887 1.17409 -2.46544
H 3.80701 1.13585 -2.60280
H 3.96018 1.03742 1.93956
H 4.10265 -0.68614 2.44380
H 1.85627 -0.45697 0.70921
H 1.90589 1.03920 1.57171
H 5.93727 5.51070 -0.50018
H 7.23660 6.54387 -0.99366
H 6.78341 6.54444 0.72721
H 9.00674 7.66464 -0.22824
H 10.39253 8.16983 0.89918
H 8.89700 7.39732 1.49919
H 13.68248 2.18867 0.68697
H 13.67972 -0.53214 -0.51885
H 14.63696 2.58095 -1.32902
H 13.03364 2.23830 -2.06928
H 14.15928 0.98803 -1.90991
H 15.22562 0.19203 1.34061
H 14.45877 -1.36679 1.57939
H 14.02814 1.25899 2.80494
H 14.14196 -0.33488 3.63647
H 12.61623 0.19761 2.79806
H 10.25576 -5.37292 1.63993
H 11.63146 -4.43175 0.97823
H 10.36035 -5.04425 -0.15375
H 5.48084 -3.46891 -0.03817
H 2.24643 -4.94546 2.39929
H 1.98777 -3.08806 2.48345
H 3.01832 -4.03389 3.58599

