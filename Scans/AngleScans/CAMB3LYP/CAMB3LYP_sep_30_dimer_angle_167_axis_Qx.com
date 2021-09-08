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
Mg 15.75805 1.41742 1.20404
C 16.08891 4.08266 3.58439
C 16.10368 -0.75404 3.81170
C 15.92170 -0.94163 -0.91150
C 16.57767 3.72044 -1.27992
N 16.05167 1.74332 3.37195
C 16.10012 2.79946 4.16980
C 16.08324 2.46228 5.68492
C 16.66220 0.96653 5.60486
C 16.21810 0.60804 4.18820
C 18.17960 0.90480 5.65974
C 14.63422 2.48592 6.28642
C 14.57229 2.47101 7.81537
H 13.51204 3.22689 8.47642
N 15.49216 -0.53564 1.43519
C 15.65417 -1.25043 2.56837
C 15.35704 -2.64588 2.20976
C 15.16686 -2.72108 0.73568
C 15.54838 -1.37678 0.35176
C 15.19950 -3.73307 3.23280
C 14.72809 -3.82455 -0.17866
O 14.86275 -3.73885 -1.45041
C 14.22236 -5.04634 0.46395
N 16.39215 1.33921 -0.77438
C 16.30942 0.19906 -1.48230
C 16.56132 0.44887 -2.95013
C 16.59195 2.03918 -3.11622
C 16.43442 2.45127 -1.63610
C 17.85110 -0.14643 -3.51448
C 15.62375 2.65630 -4.16236
C 14.19251 2.47702 -3.70981
N 16.06744 3.55339 1.10427
C 16.42870 4.31082 -0.03146
C 16.63461 5.70211 0.42057
C 16.39945 5.65261 1.80287
C 16.08997 4.32034 2.18859
C 17.03731 6.79268 -0.40078
C 16.35304 6.42242 3.02513
O 16.44231 7.61768 3.30286
C 16.27382 5.41649 4.25322
C 15.18615 5.89447 5.17573
O 13.99866 5.89425 4.80631
O 15.69766 6.06622 6.40830
C 14.65445 6.52930 7.36175
H 16.42822 -1.51706 4.51976
H 15.88613 -1.69695 -1.69999
H 16.73903 4.51686 -2.00954
H 16.73157 3.07977 6.29496
H 16.18053 0.31714 6.33126
H 18.62437 1.88820 5.62444
H 18.48160 0.32188 4.79948
H 18.50229 0.33682 6.53701
H 14.14067 1.52994 6.06905
H 14.07924 3.32333 5.91043
H 15.54500 2.69579 8.26951
H 14.35213 1.45236 8.14383
H 15.55940 -3.41314 4.21845
H 15.87844 -4.54325 2.94557
H 14.18104 -4.12471 3.27782
H 14.99096 -5.45691 1.11445
H 13.87496 -5.68327 -0.35184
H 13.37605 -4.77904 1.10654
H 15.75679 0.04788 -3.55431
H 17.58121 2.39624 -3.43825
H 17.68004 -0.82641 -4.37031
H 18.36680 -0.66275 -2.71574
H 18.59372 0.58174 -3.83618
H 15.71504 2.13344 -5.11758
H 15.80691 3.70859 -4.35244
H 13.95595 1.45997 -4.03653
H 13.54440 3.20688 -4.19075
H 14.12008 2.49717 -2.61620
H 16.42171 7.63718 -0.12086
H 16.93115 6.55515 -1.45723
H 18.08555 6.88379 -0.10386
H 17.25060 5.41436 4.73848
H 15.02156 7.45721 7.81110
H 14.47303 5.67699 8.04321
H 13.70688 6.85479 6.94867

