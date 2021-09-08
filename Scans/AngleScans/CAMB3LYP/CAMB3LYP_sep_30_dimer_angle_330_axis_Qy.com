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
Mg 15.76817 1.42421 1.20692
C 14.36636 -0.95181 3.50225
C 13.54165 3.69836 2.43389
C 16.66879 3.40397 -1.10324
C 17.05156 -1.27441 -0.58642
N 14.20689 1.27353 2.76640
C 13.80682 0.33806 3.61444
C 12.80576 0.83294 4.69265
C 12.21177 2.10045 3.90558
C 13.40289 2.40578 2.99966
C 11.01923 1.75898 3.02776
C 13.52223 1.26518 6.02000
C 12.59100 1.46524 7.21771
H 13.07627 1.09909 8.54546
N 15.55203 3.39021 1.04361
C 14.60715 4.15630 1.62808
C 14.86823 5.53175 1.17627
C 15.94288 5.49155 0.14757
C 16.08487 4.06027 -0.02952
C 14.18417 6.72767 1.77230
C 16.70634 6.56403 -0.56877
O 17.42809 6.30518 -1.59585
C 16.51019 7.93998 -0.08942
N 16.54022 1.11021 -0.69773
C 16.89724 2.14521 -1.47826
C 17.67807 1.67619 -2.68271
C 17.98068 0.12760 -2.42216
C 17.21201 -0.06078 -1.09586
C 16.97946 1.83785 -4.03260
C 19.46694 -0.31933 -2.48500
C 20.23700 0.28586 -1.33341
N 15.89277 -0.71398 1.49101
C 16.44895 -1.66136 0.60375
C 16.19672 -2.99805 1.18000
C 15.48550 -2.73469 2.36038
C 15.28948 -1.33402 2.49900
C 16.56686 -4.24058 0.59243
C 14.84674 -3.32375 3.51523
O 14.76702 -4.46085 3.97824
C 13.98370 -2.20389 4.24150
C 14.28442 -2.26779 5.71392
O 15.41940 -1.99024 6.13983
O 13.13397 -2.43283 6.39185
C 13.37922 -2.49096 7.85756
H 12.73844 4.41999 2.58545
H 17.09470 4.04913 -1.87515
H 17.50523 -2.16061 -1.03532
H 12.01073 0.13575 4.92850
H 12.02302 2.93140 4.58030
H 10.84089 0.69473 2.98657
H 11.25918 2.14074 2.04402
H 10.13635 2.31318 3.35901
H 13.90260 2.28831 5.90478
H 14.29822 0.57085 6.27762
H 11.59560 1.03971 7.04150
H 12.40724 2.53549 7.33944
H 13.32663 6.43411 2.39023
H 13.74147 7.28933 0.94258
H 14.87249 7.38039 2.31345
H 15.45619 8.19893 -0.15641
H 17.20602 8.55236 -0.66621
H 16.77773 7.98205 0.97231
H 18.61755 2.20869 -2.76588
H 17.48655 -0.51485 -3.16574
H 17.56137 2.43976 -4.75590
H 16.00770 2.28220 -3.86298
H 16.72315 0.90636 -4.53443
H 19.93537 0.05090 -3.40023
H 19.59503 -1.39654 -2.47069
H 20.48415 1.28646 -1.70057
H 21.13565 -0.29220 -1.12766
H 19.59647 0.40700 -0.45211
H 16.97100 -4.85287 1.38759
H 17.28885 -4.10172 -0.20957
H 15.59617 -4.57582 0.21741
H 12.93400 -2.40989 4.02873
H 12.94261 -3.42766 8.21764
H 12.96376 -1.55123 8.26748
H 14.40540 -2.59427 8.19031

