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
Mg 8.92557 0.80486 0.68477
C 9.98224 3.18470 3.15436
C 8.49362 -1.42221 3.23196
C 8.42447 -1.44285 -1.49791
C 10.50142 2.79324 -1.70453
N 9.23111 0.97717 2.86862
C 9.57625 1.94895 3.69981
C 9.40410 1.60119 5.20281
C 9.49481 0.00173 5.09147
C 9.01045 -0.17103 3.65334
C 10.91603 -0.52832 5.18374
C 8.01362 2.05977 5.76697
C 7.89778 2.03201 7.29261
H 7.10107 3.06487 7.94906
N 8.06132 -0.97416 0.84831
C 7.95557 -1.72810 1.96261
C 7.25417 -2.95476 1.55334
C 7.10066 -2.93569 0.07321
C 7.89190 -1.76789 -0.25885
C 6.73338 -3.96129 2.53777
C 6.37392 -3.82902 -0.88597
O 6.57195 -3.76201 -2.15053
C 5.49364 -4.84744 -0.29486
N 9.57183 0.57644 -1.27801
C 9.16506 -0.46649 -2.02284
C 9.53197 -0.27562 -3.47520
C 10.05827 1.23009 -3.59106
C 9.98522 1.63887 -2.10346
C 10.59314 -1.22904 -4.02409
C 9.36463 2.13914 -4.64233
C 7.93344 2.40253 -4.23282
N 9.88320 2.74144 0.65930
C 10.49954 3.37381 -0.44266
C 11.10976 4.62278 0.05747
C 10.82361 4.61896 1.43100
C 10.10451 3.44022 1.76710
C 11.85769 5.55223 -0.71897
C 10.97557 5.33887 2.67492
O 11.42032 6.44138 2.99179
C 10.54731 4.38091 3.86878
C 9.62969 5.15251 4.77714
O 8.51366 5.52819 4.37723
O 10.12671 5.13083 6.02729
C 9.24568 5.87383 6.96727
H 8.54200 -2.26324 3.92422
H 8.18422 -2.13288 -2.31008
H 10.92592 3.51591 -2.40468
H 10.19024 1.97415 5.84830
H 8.81144 -0.48185 5.78466
H 11.64389 0.26937 5.19048
H 11.05236 -1.15755 4.31394
H 11.01712 -1.18702 6.05105
H 7.25652 1.30868 5.50736
H 7.75785 3.03580 5.40286
H 8.87619 1.93452 7.77857
H 7.36245 1.12492 7.58352
H 7.14056 -3.78984 3.54186
H 7.13811 -4.93562 2.24326
H 5.64283 -4.01892 2.54420
H 6.07496 -5.48981 0.36231
H 4.99455 -5.32771 -1.13874
H 4.74989 -4.34488 0.33354
H 8.66408 -0.39452 -4.11217
H 11.11987 1.26985 -3.87612
H 10.24971 -1.80404 -4.90473
H 10.89634 -1.89676 -3.22863
H 11.53501 -0.76012 -4.30370
H 9.32251 1.63434 -5.61061
H 9.87045 3.08670 -4.79489
H 7.40546 1.51610 -4.59693
H 7.55951 3.30743 -4.70755
H 7.83337 2.42068 -3.14135
H 11.52403 6.53971 -0.42907
H 11.71956 5.38200 -1.78467
H 12.87195 5.30765 -0.39229
H 11.45839 4.06580 4.37890
H 9.86599 6.63244 7.45452
H 8.78650 5.10527 7.61705
H 8.45984 6.48570 6.54004

