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
Mg 13.12997 1.18439 1.00405
C 14.44504 3.38816 3.51265
C 12.38920 -0.99580 3.52053
C 12.37892 -0.96092 -1.20976
C 14.97352 2.98254 -1.34407
N 13.42656 1.29540 3.19311
C 13.88019 2.20787 4.03925
C 13.64642 1.86968 5.53610
C 13.53713 0.27268 5.40758
C 13.05348 0.17632 3.96208
C 14.87929 -0.43260 5.51018
C 12.31733 2.49371 6.08938
C 12.17928 2.46567 7.61316
H 11.51010 3.58390 8.27211
N 12.04727 -0.47349 1.13731
C 11.83341 -1.21914 2.24164
C 10.98895 -2.34387 1.81034
C 10.85812 -2.29107 0.32885
C 11.79386 -1.22869 0.01929
C 10.33334 -3.28667 2.77715
C 10.03744 -3.07652 -0.64877
O 10.25858 -3.02245 -1.91012
C 9.02877 -3.98208 -0.07952
N 13.76768 0.89598 -0.95362
C 13.24288 -0.08019 -1.71494
C 13.64951 0.07739 -3.16069
C 14.36205 1.50613 -3.25330
C 14.32172 1.90616 -1.76208
C 14.58963 -0.99632 -3.70815
C 13.80156 2.50544 -4.30207
C 12.40958 2.94248 -3.90612
N 14.32329 2.98549 1.01173
C 15.02826 3.54627 -0.07577
C 15.78390 4.70370 0.44558
C 15.48186 4.72230 1.81559
C 14.61625 3.63997 2.12990
C 16.65248 5.53944 -0.31153
C 15.70693 5.40511 3.06931
O 16.28239 6.43984 3.40384
C 15.14650 4.49680 4.24714
C 14.32134 5.36855 5.15360
O 13.26651 5.88538 4.74514
O 14.79557 5.27228 6.40907
C 14.00271 6.11076 7.34723
H 12.32274 -1.84303 4.20365
H 12.06446 -1.60726 -2.03247
H 15.49433 3.65304 -2.03099
H 14.46477 2.13455 6.19483
H 12.78960 -0.12805 6.08728
H 15.70135 0.26724 5.53442
H 14.94677 -1.06535 4.63490
H 14.88575 -1.10730 6.37104
H 11.47537 1.84627 5.81250
H 12.19077 3.49768 5.73349
H 13.13138 2.24125 8.10924
H 11.53066 1.63019 7.88752
H 10.74585 -3.17765 3.78777
H 10.61636 -4.30118 2.47624
H 9.24418 -3.20692 2.77035
H 9.51638 -4.69884 0.57695
H 8.48427 -4.38750 -0.93452
H 8.34594 -3.39630 0.54596
H 12.78183 0.07475 -3.80895
H 15.42383 1.41503 -3.52562
H 14.18813 -1.51490 -4.59919
H 14.79637 -1.70464 -2.91691
H 15.58642 -0.64669 -3.97151
H 13.70890 2.01951 -5.27646
H 14.42421 3.38342 -4.43796
H 11.77927 2.13303 -4.28637
H 12.15830 3.89182 -4.37478
H 12.29853 2.96228 -2.81574
H 16.44167 6.55811 -0.01425
H 16.50782 5.39844 -1.38062
H 17.62374 5.16618 0.02403
H 16.00419 4.06472 4.76411
H 14.70698 6.78059 7.85023
H 13.44238 5.39960 7.98287
H 13.30543 6.82068 6.91798

