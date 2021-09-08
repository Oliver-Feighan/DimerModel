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
Mg 8.92490 0.80670 0.67424
C 8.49553 -1.98918 -1.53434
C 8.80901 2.81149 -2.08293
C 9.74959 3.24235 2.53297
C 10.13214 -1.43566 3.05333
N 8.67942 0.35217 -1.47644
C 8.46310 -0.74157 -2.19149
C 8.11735 -0.47787 -3.68162
C 8.80635 0.96438 -3.83676
C 8.73280 1.42800 -2.38338
C 10.26930 0.89205 -4.24158
C 6.57015 -0.40530 -3.93244
C 6.15401 -0.45907 -5.40407
H 4.91652 -1.15161 -5.75260
N 8.75534 2.76195 0.37961
C 8.69977 3.40442 -0.80590
C 8.59707 4.83603 -0.48323
C 8.76250 4.99859 0.98678
C 9.12421 3.64651 1.36246
C 8.28451 5.88200 -1.51360
C 8.63104 6.17878 1.90115
O 9.05245 6.14338 3.11110
C 8.07954 7.40718 1.31090
N 10.00759 0.92545 2.44512
C 10.17578 2.10138 3.07505
C 10.74473 1.90203 4.45959
C 10.69740 0.32493 4.72056
C 10.16889 -0.14312 3.34688
C 12.17049 2.41027 4.67244
C 9.95808 -0.15490 5.99955
C 8.47769 0.12527 5.87570
N 9.09263 -1.34068 0.84359
C 9.65305 -2.07066 1.91464
C 9.64578 -3.49477 1.52185
C 9.09846 -3.49205 0.22998
C 8.80531 -2.15825 -0.16313
C 10.14868 -4.57504 2.30045
C 8.71187 -4.31320 -0.89472
O 8.64650 -5.52373 -1.10469
C 8.42143 -3.36478 -2.13671
C 7.11651 -3.79106 -2.75151
O 6.05143 -3.67037 -2.12091
O 7.31203 -4.06578 -4.05407
C 6.04407 -4.48258 -4.71010
H 9.01384 3.50847 -2.89604
H 9.95430 4.03517 3.25616
H 10.40103 -2.20682 3.77838
H 8.55849 -1.17786 -4.38105
H 8.21692 1.61705 -4.47548
H 10.63716 -0.12324 -4.24316
H 10.80553 1.48762 -3.51441
H 10.41860 1.38692 -5.20550
H 6.21200 0.59914 -3.67249
H 6.05888 -1.17246 -3.38408
H 6.97469 -0.78876 -6.05271
H 5.93800 0.55771 -5.74098
H 8.37992 5.48487 -2.53173
H 9.06893 6.64343 -1.44445
H 7.31502 6.35753 -1.35044
H 8.70273 7.71785 0.47569
H 7.97967 8.11045 2.13995
H 7.08938 7.18323 0.89855
H 10.13500 2.39979 5.20361
H 11.70582 -0.10027 4.83077
H 12.25408 3.14309 5.49710
H 12.52153 2.84086 3.74405
H 12.91265 1.63709 4.86369
H 10.30781 0.40371 6.87121
H 10.10343 -1.20865 6.21280
H 8.39871 1.17356 6.17891
H 7.90847 -0.52174 6.54009
H 8.15052 0.05857 4.83159
H 9.42490 -5.37574 2.22620
H 10.30984 -4.27840 3.33477
H 11.08906 -4.77063 1.77824
H 9.25531 -3.47052 -2.83171
H 6.22747 -5.45937 -5.16827
H 5.77105 -3.65173 -5.38753
H 5.19811 -4.70461 -4.07018

