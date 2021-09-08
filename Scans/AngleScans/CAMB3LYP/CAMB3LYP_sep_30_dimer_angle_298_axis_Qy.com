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
Mg 15.76258 1.42371 1.20661
C 13.52699 -1.15065 2.32636
C 13.06786 3.51247 1.10564
C 17.61251 3.58052 -0.20526
C 17.99605 -1.09771 0.31235
N 13.62535 1.10555 1.67912
C 12.90279 0.10247 2.15459
C 11.44797 0.48345 2.53922
C 11.27674 1.76750 1.59013
C 12.74193 2.18273 1.47399
C 10.76580 1.42493 0.20055
C 11.31125 0.86690 4.05457
C 9.87266 0.94989 4.56966
H 9.59853 0.52650 5.94008
N 15.52605 3.38379 1.00483
C 14.36276 4.06790 1.01049
C 14.72488 5.47902 0.80610
C 16.18246 5.55180 0.51476
C 16.49981 4.14217 0.40385
C 13.74411 6.60316 0.97293
C 17.13086 6.70041 0.34910
O 18.30609 6.53862 -0.13616
C 16.61040 8.03392 0.68432
N 17.45294 1.26280 0.00697
C 18.09599 2.35900 -0.43207
C 19.43088 2.00123 -1.04072
C 19.65821 0.45545 -0.69947
C 18.31586 0.15088 0.00106
C 19.55189 2.21214 -2.54978
C 20.97585 0.08300 0.03385
C 20.96569 0.65210 1.43435
N 15.86990 -0.71967 1.45652
C 16.88060 -1.58420 0.98203
C 16.45668 -2.96385 1.29712
C 15.20810 -2.80642 1.91742
C 14.86799 -1.42749 1.96606
C 17.17187 -4.14962 0.96761
C 14.09564 -3.49397 2.53256
O 13.86312 -4.65966 2.85008
C 12.90003 -2.46196 2.71109
C 12.37139 -2.60142 4.11237
O 13.07985 -2.29772 5.08832
O 11.05213 -2.86099 4.06183
C 10.47994 -2.99688 5.42783
H 12.25835 4.18526 0.82118
H 18.33709 4.29071 -0.61034
H 18.68185 -1.93276 0.15405
H 10.70260 -0.26316 2.29272
H 10.69744 2.54605 2.07964
H 10.71430 0.35801 0.04191
H 11.46588 1.87674 -0.48996
H 9.80527 1.91611 0.02040
H 11.61943 1.91177 4.18873
H 11.87729 0.19510 4.67014
H 9.15904 0.48993 3.87530
H 9.57572 2.00099 4.60191
H 12.71288 6.23278 1.02570
H 13.77396 7.19310 0.05047
H 13.98775 7.25306 1.81610
H 15.73975 8.24702 0.06869
H 17.46047 8.71212 0.58784
H 16.26556 8.02360 1.72432
H 20.22824 2.58106 -0.59226
H 19.68551 -0.16282 -1.60877
H 20.38501 2.88352 -2.83125
H 18.61074 2.59949 -2.91685
H 19.67120 1.30204 -3.13518
H 21.83262 0.52980 -0.47651
H 21.15374 -0.98593 0.08621
H 21.29786 1.68367 1.28395
H 21.65437 0.10513 2.07512
H 19.94667 0.68927 1.83668
H 17.13165 -4.78974 1.83890
H 18.19840 -3.92858 0.68277
H 16.57876 -4.50629 0.12128
H 12.14438 -2.70335 1.96255
H 9.98736 -3.97314 5.47222
H 9.84329 -2.10429 5.57520
H 11.17398 -3.07230 6.25673

