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
Mg 10.50520 0.94847 0.78975
C 11.52485 -1.72451 -1.37694
C 11.73591 3.11006 -1.54423
C 10.03136 3.30206 2.86424
C 10.41281 -1.37602 3.38473
N 11.48005 0.61497 -1.16746
C 11.75865 -0.44203 -1.91546
C 12.24442 -0.10438 -3.35058
C 12.80367 1.37522 -3.07327
C 11.93194 1.74542 -1.87506
C 14.25736 1.39603 -2.63091
C 11.06997 -0.08890 -4.39099
C 11.50946 -0.07220 -5.85671
H 10.70322 -0.79917 -6.83342
N 10.37883 2.90800 0.50198
C 10.91899 3.61828 -0.51046
C 10.55694 5.02119 -0.25569
C 9.89931 5.10136 1.07691
C 10.10075 3.74723 1.55215
C 10.76868 6.11232 -1.26466
C 9.21507 6.21617 1.80848
O 8.92624 6.12676 3.05404
C 8.97731 7.45121 1.04713
N 10.46261 1.00939 2.86745
C 10.18301 2.15129 3.52000
C 9.93701 1.89465 4.98766
C 9.87142 0.30408 5.14058
C 10.19379 -0.10348 3.68615
C 10.98758 2.45492 5.94619
C 8.60006 -0.28681 5.80928
C 7.39910 -0.06899 4.91729
N 10.71079 -1.19508 0.96584
C 10.66331 -1.96206 2.15056
C 10.96966 -3.35836 1.77781
C 11.19855 -3.30241 0.39471
C 11.06549 -1.96225 -0.05890
C 11.05520 -4.45945 2.67581
C 11.53292 -4.07059 -0.78291
O 11.67727 -5.26779 -1.02706
C 11.88336 -3.06278 -1.96100
C 11.14338 -3.51119 -3.19139
O 9.90072 -3.47899 -3.22845
O 12.02367 -3.69655 -4.19194
C 11.33616 -4.13127 -5.43700
H 12.29255 3.86412 -2.10139
H 9.76033 4.05800 3.60490
H 10.30757 -2.17657 4.12007
H 13.04006 -0.73907 -3.72193
H 12.60135 2.03741 -3.91108
H 14.64124 0.40100 -2.46153
H 14.27766 1.97052 -1.71407
H 14.86236 1.95519 -3.35039
H 10.55705 0.88001 -4.33758
H 10.40171 -0.91109 -4.22347
H 12.57097 -0.32307 -5.97163
H 11.43420 0.95205 -6.22990
H 11.42154 5.78289 -2.08238
H 11.33756 6.90388 -0.76499
H 9.83065 6.53127 -1.63513
H 9.92598 7.84097 0.68577
H 8.39963 8.09719 1.71111
H 8.37972 7.20687 0.16181
H 8.99012 2.31713 5.30072
H 10.69257 -0.07954 5.76378
H 10.56465 3.13919 6.70578
H 11.74814 2.95723 5.36330
H 11.56627 1.70698 6.48555
H 8.38875 0.23331 6.74682
H 8.68453 -1.34368 6.03928
H 7.09505 0.95404 5.15825
H 6.61139 -0.78115 5.15471
H 7.68608 -0.08708 3.85948
H 10.54294 -5.28703 2.20336
H 10.61704 -4.21925 3.64229
H 12.14032 -4.57876 2.73530
H 12.96462 -3.08693 -2.10201
H 11.80579 -5.06869 -5.75065
H 11.40819 -3.27431 -6.13283
H 10.29779 -4.43114 -5.35759

