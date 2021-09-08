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
Mg 31.54176 2.84280 2.39191
C 33.87699 0.70458 4.08123
C 32.97831 5.46232 4.03783
C 29.21880 4.69145 1.27218
C 29.59879 0.01293 1.78990
N 33.20600 2.94359 3.84532
C 34.05896 2.07405 4.36562
C 35.21773 2.71542 5.17522
C 34.50499 4.10283 5.55730
C 33.52271 4.20204 4.39198
C 33.70792 4.03414 6.84931
C 36.49948 2.94977 4.30104
C 37.76560 3.28416 5.09281
H 39.04606 2.78142 4.60278
N 31.43979 4.81771 2.22486
C 32.10843 5.73200 2.95852
C 31.67073 7.04534 2.46079
C 30.56648 6.83329 1.48592
C 30.34670 5.40944 1.64218
C 32.34197 8.32933 2.85352
C 29.82178 7.76106 0.57440
O 28.73958 7.40382 -0.01218
C 30.34559 9.13018 0.46204
N 29.57962 2.45396 1.82520
C 28.79939 3.42699 1.32283
C 27.53122 2.85876 0.73185
C 27.73354 1.27237 0.72525
C 29.10717 1.18953 1.42643
C 26.23938 3.20372 1.47246
C 27.55328 0.54556 -0.63560
C 28.65989 0.94537 -1.58479
N 31.76681 0.71266 2.67458
C 30.82014 -0.29201 2.37701
C 31.38217 -1.57368 2.85016
C 32.61657 -1.21972 3.41538
C 32.80118 0.18613 3.32058
C 30.74096 -2.84216 2.77311
C 33.80065 -1.71434 4.08020
O 34.24190 -2.83005 4.35263
C 34.61292 -0.47462 4.65439
C 36.05792 -0.65034 4.27542
O 36.40706 -0.61638 3.08229
O 36.81317 -0.61166 5.38829
C 38.25562 -0.77062 5.06287
H 33.20403 6.32125 4.67043
H 28.43321 5.26863 0.77905
H 29.09803 -0.92783 1.55142
H 35.49425 2.17928 6.07514
H 35.21053 4.92925 5.53267
H 33.65523 3.02769 7.23717
H 32.71849 4.39744 6.60423
H 34.11472 4.73798 7.58106
H 36.38070 3.88210 3.73411
H 36.68436 2.10787 3.66256
H 37.65196 3.06865 6.16220
H 37.92488 4.36445 5.05256
H 33.01328 8.18690 3.70935
H 31.55953 8.99698 3.23021
H 32.84690 8.81419 2.01529
H 30.36084 9.59287 1.44599
H 29.73432 9.61519 -0.30149
H 31.38592 9.07991 0.12174
H 27.39273 3.19959 -0.28685
H 27.01297 0.76698 1.38495
H 25.49000 3.70618 0.83212
H 26.48890 3.82327 2.32360
H 25.73605 2.35960 1.94080
H 26.61538 0.85029 -1.10639
H 27.53350 -0.53565 -0.54853
H 28.29899 1.89099 -2.00028
H 28.78686 0.19530 -2.36281
H 29.58770 1.15712 -1.04074
H 31.49056 -3.55060 2.44644
H 29.89234 -2.81799 2.09266
H 30.42938 -2.96671 3.81362
H 34.47181 -0.46303 5.73583
H 38.62462 -1.61625 5.65146
H 38.71598 0.21656 5.25597
H 38.51095 -1.08513 4.05771

