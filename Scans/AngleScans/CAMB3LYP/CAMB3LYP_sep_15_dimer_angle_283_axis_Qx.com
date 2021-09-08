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
Mg 7.88828 0.70469 0.60663
C 5.26434 0.35505 3.02981
C 9.60488 2.49936 2.94480
C 10.01299 1.41265 -1.64098
C 5.54346 -0.09838 -1.84339
N 7.39390 1.28446 2.68310
C 6.39493 1.04239 3.51845
C 6.65441 1.49611 4.98008
C 7.75261 2.63086 4.68824
C 8.31613 2.09349 3.37433
C 7.14797 3.99982 4.42423
C 7.24141 0.34360 5.86842
C 7.22982 0.62197 7.37312
H 6.98686 -0.49224 8.28525
N 9.75420 1.36191 0.76458
C 10.28902 2.05362 1.79259
C 11.68647 2.31746 1.41620
C 11.88671 1.86851 0.01154
C 10.52388 1.54225 -0.35773
C 12.69895 2.86581 2.37934
C 13.10160 1.76382 -0.85982
O 13.00916 1.57361 -2.12399
C 14.39924 1.98361 -0.20460
N 7.75071 0.87226 -1.46105
C 8.83714 1.13963 -2.20678
C 8.55273 0.93973 -3.67632
C 7.11328 0.24630 -3.74490
C 6.76038 0.24655 -2.24128
C 8.54826 2.20670 -4.53139
C 7.01323 -1.09175 -4.52739
C 7.78026 -2.17307 -3.80090
N 5.83433 0.03489 0.57946
C 5.03420 -0.21685 -0.55654
C 3.67689 -0.53895 -0.07034
C 3.77659 -0.42888 1.32481
C 5.09726 -0.04424 1.68170
C 2.54816 -0.83072 -0.88702
C 3.06215 -0.54026 2.57616
O 1.93787 -0.91743 2.90397
C 3.96053 0.08553 3.72841
C 3.97080 -0.87727 4.88403
O 4.49875 -1.99736 4.76882
O 3.55241 -0.23964 5.99252
C 3.55575 -1.15043 7.16822
H 10.12702 3.26365 3.52111
H 10.73522 1.53847 -2.45094
H 4.78116 -0.45555 -2.53910
H 5.79647 1.92892 5.48047
H 8.52161 2.64313 5.45626
H 6.07081 3.95883 4.35961
H 7.57095 4.33211 3.48524
H 7.48958 4.71459 5.17838
H 8.32363 0.27608 5.69817
H 6.74101 -0.58394 5.66867
H 6.59000 1.47509 7.62956
H 8.23133 0.93767 7.67534
H 12.22130 3.25398 3.28740
H 13.14355 3.74880 1.90762
H 13.49191 2.15095 2.60925
H 14.41333 2.97324 0.24583
H 15.15109 1.77915 -0.96932
H 14.50396 1.26346 0.61457
H 9.28333 0.27796 -4.12509
H 6.37352 0.89133 -4.24141
H 9.26373 2.16689 -5.37436
H 8.76227 3.05174 -3.89057
H 7.58138 2.47571 -4.95327
H 7.47679 -0.99037 -5.51177
H 5.99280 -1.42485 -4.68473
H 8.80967 -2.01227 -4.13489
H 7.42066 -3.15836 -4.09063
H 7.75644 -2.01189 -2.71688
H 2.04542 -1.67450 -0.43338
H 2.84392 -1.04697 -1.91148
H 2.00125 0.11264 -0.80810
H 3.52219 1.04657 4.00013
H 2.54580 -1.13220 7.58935
H 4.37871 -0.79706 7.81753
H 3.68783 -2.20970 6.98124

