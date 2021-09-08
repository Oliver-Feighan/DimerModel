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
Mg 7.89060 0.71125 0.60619
C 6.60619 -1.27876 3.30234
C 7.69365 3.40040 2.69528
C 8.58977 2.56328 -1.87343
C 6.91860 -1.84218 -1.55734
N 7.22105 0.90852 2.70501
C 6.81241 0.07305 3.64806
C 6.69048 0.70664 5.05989
C 6.49607 2.23879 4.61991
C 7.21024 2.19808 3.27047
C 5.04372 2.61652 4.38047
C 7.99191 0.51537 5.91530
C 7.83575 0.83427 7.40381
H 8.58313 0.03328 8.36938
N 8.57112 2.57455 0.54645
C 8.42051 3.52654 1.49108
C 9.08485 4.72750 0.96146
C 9.49668 4.45517 -0.44244
C 8.86964 3.16621 -0.65577
C 9.34461 5.94833 1.79545
C 10.30841 5.23274 -1.43375
O 10.34117 4.91416 -2.67482
C 10.98905 6.43296 -0.92612
N 7.58291 0.50410 -1.43947
C 8.03133 1.43168 -2.30345
C 7.94098 0.93911 -3.72809
C 7.56241 -0.61116 -3.62348
C 7.40502 -0.73041 -2.09168
C 6.92140 1.65295 -4.61537
C 8.49858 -1.61485 -4.35069
C 9.85186 -1.63921 -3.67731
N 7.10653 -1.29105 0.81616
C 6.74440 -2.17840 -0.22093
C 6.15534 -3.37448 0.41533
C 6.19488 -3.08797 1.78823
C 6.74956 -1.79556 1.99196
C 5.63016 -4.50705 -0.26846
C 5.88373 -3.58058 3.11081
O 5.47778 -4.64881 3.56668
C 6.01972 -2.37912 4.14251
C 6.82116 -2.86545 5.31882
O 8.01593 -3.18234 5.18200
O 6.11235 -2.67056 6.44575
C 6.87010 -3.12722 7.64121
H 7.45903 4.34326 3.19011
H 8.91415 3.11540 -2.75855
H 6.68130 -2.72309 -2.15764
H 5.83486 0.37349 5.63496
H 7.00733 2.91344 5.30169
H 4.39053 1.75809 4.43124
H 5.01338 3.05711 3.39254
H 4.74088 3.40744 5.07262
H 8.72180 1.28527 5.63377
H 8.38305 -0.47645 5.79733
H 6.78268 0.91032 7.70080
H 8.23922 1.83274 7.58898
H 8.78222 5.91948 2.73688
H 8.92298 6.80162 1.25323
H 10.40923 6.12626 1.96156
H 10.25306 7.11648 -0.50947
H 11.58939 6.80295 -1.75953
H 11.64795 6.13906 -0.10145
H 8.89518 1.03464 -4.23155
H 6.57369 -0.81937 -4.05812
H 7.36803 2.09156 -5.52763
H 6.43212 2.41722 -4.02628
H 6.08274 1.03970 -4.94061
H 8.67037 -1.29411 -5.38112
H 8.10296 -2.62430 -4.38933
H 10.36465 -0.78152 -4.12283
H 10.37280 -2.56873 -3.89786
H 9.75998 -1.44630 -2.60214
H 5.98387 -5.38188 0.26075
H 5.93940 -4.52038 -1.31152
H 4.55836 -4.32040 -0.16044
H 5.01150 -2.07862 4.43017
H 6.23510 -3.84773 8.16594
H 7.14670 -2.20687 8.18911
H 7.76406 -3.71630 7.47306

