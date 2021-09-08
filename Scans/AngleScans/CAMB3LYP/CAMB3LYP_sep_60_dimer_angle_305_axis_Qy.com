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
Mg 31.53007 2.84325 2.40462
C 29.41873 0.31227 3.82428
C 28.87190 4.96994 2.61900
C 33.21664 4.95961 0.74814
C 33.60005 0.28135 1.26554
N 29.46499 2.56065 3.14479
C 28.79457 1.57267 3.71819
C 27.40470 1.97836 4.27786
C 27.13276 3.25582 3.34352
C 28.57717 3.64859 3.04004
C 26.44700 2.90773 2.03283
C 27.46430 2.37792 5.79405
C 26.10294 2.48664 6.48444
H 25.99736 2.08013 7.88303
N 31.29610 4.80459 2.21231
C 30.15188 5.50556 2.35613
C 30.50416 6.90929 2.09231
C 31.91450 6.95814 1.61978
C 32.19671 5.54306 1.48563
C 29.56707 8.04909 2.36825
C 32.84975 8.09123 1.32380
O 33.95256 7.90783 0.69690
C 32.39320 9.43524 1.70678
N 33.05425 2.64655 1.00443
C 33.65160 3.72912 0.47611
C 34.89471 3.34630 -0.29105
C 35.14254 1.80063 0.03617
C 33.89479 1.52219 0.90277
C 34.82818 3.54128 -1.80562
C 36.53674 1.41593 0.60265
C 36.70994 1.99821 1.98697
N 31.63947 0.70097 2.66290
C 32.57110 -0.18258 2.07511
C 32.17182 -1.55291 2.45616
C 31.01314 -1.37153 3.22627
C 30.70013 0.01259 3.30184
C 32.82422 -2.75199 2.05280
C 29.97764 -2.03704 3.98360
O 29.77137 -3.19621 4.34068
C 28.82767 -0.98612 4.29909
C 28.47725 -1.10476 5.75705
O 29.30653 -0.80227 6.63306
O 27.15875 -1.34557 5.87515
C 26.76073 -1.46034 7.30352
H 28.04209 5.65175 2.43079
H 33.89401 5.65537 0.24757
H 34.24943 -0.56503 1.03181
H 26.62445 1.24039 4.13504
H 26.62984 4.04726 3.89310
H 26.36186 1.84023 1.89376
H 27.06083 3.34283 1.25505
H 25.47805 3.41112 1.96905
H 27.80070 3.41941 5.87691
H 28.09415 1.70375 6.34124
H 25.30180 2.03063 5.89015
H 25.82635 3.54220 6.54197
H 28.54579 7.69426 2.55396
H 29.48875 8.62988 1.44282
H 29.92318 8.70325 3.16696
H 31.45508 9.65522 1.20280
H 33.23336 10.10008 1.49699
H 32.18148 9.43968 2.78185
H 35.74965 3.91864 0.04747
H 35.04734 1.17354 -0.86246
H 35.62823 4.19782 -2.19673
H 33.85363 3.93884 -2.05611
H 34.86102 2.62411 -2.39121
H 37.32854 1.84543 -0.01601
H 36.70560 0.34507 0.64417
H 37.03426 3.02338 1.78467
H 37.46627 1.44730 2.54241
H 35.75004 2.05396 2.51341
H 32.88517 -3.38327 2.92931
H 33.80974 -2.54858 1.63910
H 32.12493 -3.10789 1.29149
H 27.98094 -1.22349 3.65390
H 26.26472 -2.42887 7.42014
H 26.15950 -0.55725 7.51960
H 27.55224 -1.53808 8.03967

