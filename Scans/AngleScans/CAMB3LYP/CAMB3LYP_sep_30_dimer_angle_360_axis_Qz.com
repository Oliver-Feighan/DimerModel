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
Mg 15.77261 1.42431 1.20460
C 15.55177 -0.79841 4.01348
C 14.51592 3.86588 3.22798
C 15.53076 3.28315 -1.35541
C 15.91257 -1.39534 -0.83885
N 15.15336 1.41301 3.32796
C 15.17235 0.52683 4.31217
C 14.83870 1.11170 5.71080
C 13.99404 2.39381 5.24004
C 14.61174 2.59962 3.85846
C 12.51238 2.10640 5.06279
C 16.12490 1.53139 6.50551
C 15.89726 1.82075 7.99084
H 16.94687 1.45862 8.93934
N 15.58981 3.39525 1.06125
C 15.07868 4.23448 1.98630
C 15.14845 5.57857 1.39240
C 15.59005 5.44392 -0.02237
C 15.56679 4.00245 -0.16972
C 14.88963 6.83037 2.17944
C 15.95841 6.44708 -1.07317
O 16.08266 6.11674 -2.30535
C 16.07782 7.84452 -0.63227
N 15.51480 1.01330 -0.81702
C 15.49505 2.00324 -1.72672
C 15.57624 1.45493 -3.13141
C 15.89969 -0.10257 -2.96662
C 15.85950 -0.20778 -1.42626
C 14.32049 1.62603 -3.98593
C 17.15000 -0.64443 -3.71213
C 18.40570 -0.06009 -3.10620
N 15.92620 -0.70950 1.50679
C 15.94332 -1.71267 0.51309
C 15.94305 -3.01547 1.20961
C 15.90222 -2.67737 2.57068
C 15.85831 -1.26400 2.71194
C 15.92935 -4.29362 0.58330
C 15.87555 -3.19495 3.91979
O 15.98013 -4.31248 4.42370
C 15.51977 -2.00434 4.91088
C 16.49088 -2.04951 6.05867
O 17.70125 -1.83378 5.87139
O 15.80476 -2.12378 7.21370
C 16.72436 -2.15985 8.38197
H 13.91777 4.64070 3.70835
H 15.55856 3.88000 -2.27000
H 16.05431 -2.31975 -1.40281
H 14.22700 0.47265 6.33627
H 14.19065 3.25206 5.87731
H 12.29044 1.05491 5.16887
H 12.26386 2.44686 4.06611
H 11.92406 2.72383 5.74762
H 16.44626 2.52502 6.16762
H 16.89796 0.79601 6.39496
H 14.92312 1.45491 8.33748
H 15.84162 2.90327 8.12894
H 14.42518 6.60775 3.14801
H 14.12616 7.39759 1.63608
H 15.78119 7.45170 2.28778
H 15.13478 8.16807 -0.19816
H 16.43449 8.39653 -1.50402
H 16.82623 7.89667 0.16634
H 16.38088 1.92441 -3.68390
H 15.08069 -0.73117 -3.34601
H 14.50633 2.17108 -4.93063
H 13.57164 2.13537 -3.39404
H 13.81371 0.70016 -4.25265
H 17.13389 -0.32822 -4.75803
H 17.22231 -1.72688 -3.70403
H 18.48793 0.91308 -3.59942
H 19.26603 -0.68858 -3.32727
H 18.27616 0.12399 -2.03336
H 16.64022 -4.90979 1.11762
H 16.17968 -4.22144 -0.47306
H 14.88475 -4.57604 0.73924
H 14.48999 -2.14866 5.23997
H 16.47577 -3.05742 8.95665
H 16.59954 -1.18545 8.89061
H 17.77811 -2.31955 8.18561

