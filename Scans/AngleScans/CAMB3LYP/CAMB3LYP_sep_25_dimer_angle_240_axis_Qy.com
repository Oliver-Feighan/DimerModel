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
Mg 13.12854 1.18583 0.99977
C 11.34276 -1.62435 -0.33922
C 11.51008 3.11896 -1.29767
C 15.00926 3.59584 1.84960
C 15.39275 -1.08219 2.36914
N 11.64110 0.69522 -0.56198
C 10.99124 -0.39622 -0.93719
C 9.84760 -0.14545 -1.95639
C 10.38331 1.23183 -2.58487
C 11.20313 1.73859 -1.40006
C 11.32118 1.03087 -3.76357
C 8.45555 0.04502 -1.25786
C 7.25032 -0.01435 -2.19895
H 6.01262 -0.61384 -1.70783
N 12.91502 3.13981 0.72598
C 12.20363 3.75272 -0.24332
C 12.38176 5.19607 -0.02153
C 13.38922 5.38402 1.05756
C 13.83532 4.01769 1.24258
C 11.57446 6.23591 -0.74292
C 13.88027 6.59500 1.79144
O 14.93118 6.55809 2.52427
C 13.14857 7.84767 1.55282
N 15.05089 1.26415 1.78768
C 15.61585 2.43910 2.11666
C 16.88050 2.23160 2.91541
C 16.91798 0.67101 3.26289
C 15.65930 0.21056 2.49536
C 18.18098 2.62863 2.21720
C 17.05138 0.28680 4.76198
C 15.79829 0.68210 5.50957
N 13.25728 -0.96278 1.18673
C 14.30410 -1.70669 1.77397
C 13.99634 -3.13556 1.55944
C 12.79399 -3.12325 0.83645
C 12.39201 -1.78114 0.59845
C 14.80714 -4.23168 1.96839
C 11.77898 -3.94032 0.21127
O 11.54263 -5.14680 0.16365
C 10.86003 -3.00511 -0.68723
C 9.42427 -3.34098 -0.39008
O 8.94270 -3.11791 0.73461
O 8.80103 -3.66533 -1.53757
C 7.37126 -3.99585 -1.29601
H 11.23080 3.77513 -2.12248
H 15.63965 4.38861 2.25897
H 15.99856 -1.85286 2.85052
H 9.75671 -0.89721 -2.73122
H 9.56396 1.91262 -2.80063
H 11.56669 -0.01047 -3.91011
H 12.21155 1.60047 -3.53150
H 10.89826 1.48612 -4.66370
H 8.36955 1.08168 -0.90760
H 8.32825 -0.66343 -0.46255
H 7.51380 -0.42642 -3.18064
H 6.92796 1.00717 -2.41482
H 11.03197 5.80528 -1.59361
H 12.28556 6.93313 -1.19908
H 10.91219 6.79232 -0.07624
H 13.17455 8.08458 0.49185
H 13.59122 8.57863 2.23226
H 12.09596 7.69340 1.81515
H 16.85166 2.79667 3.83897
H 17.77518 0.16892 2.79073
H 18.77048 3.37412 2.78361
H 17.93852 3.00450 1.23208
H 18.85395 1.80347 1.99045
H 17.87461 0.83848 5.22252
H 17.24207 -0.76916 4.92191
H 15.96511 1.74108 5.72824
H 15.69855 0.10115 6.42422
H 14.91603 0.61412 4.86254
H 14.14000 -4.97306 2.38764
H 15.56112 -3.92145 2.68885
H 15.24828 -4.51642 1.00944
H 11.11793 -3.19635 -1.72958
H 7.20092 -4.99616 -1.70593
H 6.79312 -3.16410 -1.74065
H 7.05494 -4.13169 -0.26843

