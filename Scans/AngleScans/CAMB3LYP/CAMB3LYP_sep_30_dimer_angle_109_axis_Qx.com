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
Mg 15.75948 1.42515 1.20341
C 17.88211 2.49963 3.89023
C 13.80477 -0.09011 3.55239
C 14.12396 0.08786 -1.16390
C 18.42778 2.02063 -0.95791
N 15.92561 1.29308 3.40505
C 16.74052 1.79264 4.32195
C 16.26486 1.58802 5.78531
C 15.32694 0.30481 5.55608
C 14.96247 0.52727 4.08968
C 16.07357 -1.01405 5.66790
C 15.44291 2.81154 6.32327
C 15.21202 2.81603 7.83602
H 15.20340 4.09815 8.53499
N 13.95191 0.60591 1.19371
C 13.30071 0.05908 2.24166
C 12.01591 -0.42183 1.71046
C 12.03079 -0.26222 0.23092
C 13.40755 0.13950 0.02296
C 10.89613 -0.89322 2.59211
C 10.98322 -0.45346 -0.82369
O 11.28088 -0.48861 -2.06989
C 9.61194 -0.69168 -0.35004
N 16.27035 0.89854 -0.74185
C 15.35587 0.38099 -1.58077
C 15.87720 0.33896 -2.99740
C 17.24773 1.16282 -2.97545
C 17.33027 1.47651 -1.46542
C 16.13058 -1.05463 -3.57193
C 17.37869 2.33794 -3.98273
C 16.41402 3.44244 -3.61548
N 17.72764 2.30137 1.36665
C 18.69253 2.42793 0.34340
C 19.91411 2.98141 0.96294
C 19.58008 3.11800 2.31877
C 18.25142 2.66163 2.53298
C 21.14232 3.24177 0.29225
C 20.05297 3.53454 3.61937
O 21.06935 4.08711 4.03814
C 19.01803 3.03468 4.71729
C 18.72999 4.18543 5.64213
O 18.14459 5.20029 5.22501
O 18.99591 3.81144 6.90695
C 18.71522 4.91583 7.86265
H 13.25093 -0.78904 4.17981
H 13.56711 -0.26293 -2.03588
H 19.27005 2.32653 -1.58215
H 17.05284 1.35146 6.49022
H 14.43846 0.34833 6.18064
H 17.13882 -0.86685 5.76647
H 15.84916 -1.55709 4.75916
H 15.66183 -1.61214 6.48587
H 14.40543 2.72684 5.97511
H 15.89659 3.73641 6.02423
H 15.86157 2.10016 8.35425
H 14.20079 2.45227 8.03379
H 11.23588 -1.05362 3.62276
H 10.61154 -1.89119 2.24150
H 10.02176 -0.24040 2.54778
H 9.59638 -1.57765 0.28024
H 8.99229 -0.71482 -1.24858
H 9.30920 0.15007 0.28274
H 15.18725 0.82270 -3.67790
H 18.11122 0.52359 -3.21127
H 15.57321 -1.24891 -4.50771
H 15.87412 -1.78661 -2.81763
H 17.17449 -1.28779 -3.77464
H 17.10440 2.00770 -4.98765
H 18.38172 2.74728 -4.04084
H 15.47493 3.11058 -4.06834
H 16.74080 4.39171 -4.03517
H 16.25979 3.48584 -2.53106
H 21.49021 4.20455 0.64229
H 21.01492 3.23303 -0.78824
H 21.73889 2.39502 0.64210
H 19.47560 2.19392 5.24010
H 19.63390 5.08657 8.43234
H 17.82125 4.59846 8.43173
H 18.53565 5.90186 7.45040

