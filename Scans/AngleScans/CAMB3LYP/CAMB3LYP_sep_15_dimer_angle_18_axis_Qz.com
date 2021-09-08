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
Mg 7.88959 0.71240 0.60680
C 7.80024 -2.27591 2.59202
C 6.65113 2.37154 3.31753
C 7.53139 3.26486 -1.24362
C 8.01868 -1.33112 -2.21332
N 7.33816 0.02704 2.63620
C 7.40521 -1.11938 3.29617
C 7.10487 -1.00577 4.81487
C 6.22139 0.33494 4.78563
C 6.79102 0.97551 3.52157
C 4.74063 0.07518 4.56381
C 8.40731 -0.81689 5.66925
C 8.22129 -1.00802 7.17611
H 9.30712 -1.61647 7.93982
N 7.66472 2.62498 1.08642
C 7.16717 3.12199 2.23825
C 7.19242 4.58530 2.08912
C 7.59153 4.90762 0.69214
C 7.59117 3.58251 0.10549
C 6.93476 5.52423 3.23168
C 7.90723 6.19668 -0.00418
O 7.99872 6.26760 -1.28067
C 8.01380 7.39176 0.84547
N 7.57585 0.94033 -1.43574
C 7.50841 2.16229 -1.99269
C 7.55558 2.07824 -3.49971
C 7.91378 0.55592 -3.83409
C 7.92437 -0.02203 -2.40187
C 6.27038 2.47025 -4.22839
C 9.15000 0.30665 -4.74090
C 10.41285 0.70960 -4.01411
N 8.09341 -1.40497 0.22831
C 8.09824 -2.05037 -1.02768
C 8.14490 -3.50427 -0.76990
C 8.14071 -3.60538 0.62947
C 8.07429 -2.30697 1.20325
C 8.13582 -4.52562 -1.76130
C 8.16663 -4.51569 1.75168
O 8.30843 -5.73093 1.88136
C 7.81970 -3.70080 3.07153
C 8.82733 -4.07187 4.12486
O 10.02682 -3.77499 3.98456
O 8.17964 -4.51920 5.21609
C 9.13626 -4.88940 6.29295
H 6.05378 2.94255 4.02887
H 7.51883 4.11605 -1.92831
H 8.16014 -2.03121 -3.03948
H 6.52558 -1.82389 5.22581
H 6.42163 0.95891 5.65277
H 4.54229 -0.96327 4.34379
H 4.45424 0.70035 3.72824
H 4.16258 0.43358 5.42040
H 8.69878 0.24111 5.64863
H 9.19040 -1.46000 5.31738
H 7.26578 -1.49030 7.41558
H 8.14938 -0.02342 7.64453
H 6.50553 4.99982 4.09433
H 6.14379 6.21020 2.90971
H 7.81727 6.10619 3.50586
H 7.07896 7.53853 1.38119
H 8.33210 8.19626 0.17963
H 8.78597 7.21505 1.60253
H 8.33322 2.71800 -3.89856
H 7.09535 0.05298 -4.36987
H 6.41579 3.28594 -4.96159
H 5.53104 2.75020 -3.48978
H 5.77320 1.65863 -4.75686
H 9.09477 0.93051 -5.63639
H 9.24320 -0.72268 -5.07075
H 10.46083 1.78955 -4.18292
H 11.27761 0.20474 -4.43996
H 10.31383 0.54884 -2.93429
H 8.87488 -5.25678 -1.46184
H 8.35118 -4.12301 -2.74888
H 7.10226 -4.87154 -1.67549
H 6.80379 -3.96865 3.36438
H 8.92317 -5.92736 6.56666
H 9.00900 -4.12421 7.08164
H 10.18613 -4.95093 6.03134

