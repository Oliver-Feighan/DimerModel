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
Mg 10.51024 0.95090 0.80796
C 8.65426 -1.51337 2.64147
C 7.98114 3.13881 1.47973
C 11.89661 3.00697 -1.17142
C 12.27979 -1.67135 -0.65431
N 8.61566 0.72416 1.92663
C 8.04577 -0.24094 2.63243
C 6.79371 0.20297 3.43547
C 6.37603 1.47359 2.54669
C 7.74516 1.83008 1.97128
C 5.44991 1.12338 1.39387
C 7.14447 0.62129 4.90652
C 5.93921 0.76944 5.83772
H 6.09009 0.38406 7.23814
N 10.28230 2.91421 0.62904
C 9.19917 3.64230 0.97254
C 9.52275 5.03421 0.62319
C 10.82005 5.04539 -0.10593
C 11.04467 3.62273 -0.26612
C 8.67630 6.19813 1.05007
C 11.70488 6.15336 -0.59138
O 12.66671 5.93716 -1.41051
C 11.35438 7.51218 -0.15292
N 11.74051 0.70171 -0.84935
C 12.24894 1.76359 -1.49883
C 13.31837 1.34304 -2.47854
C 13.59328 -0.20326 -2.17690
C 12.52515 -0.44231 -1.08724
C 12.97261 1.51924 -3.95698
C 15.06130 -0.61126 -1.87506
C 15.50241 -0.01457 -0.55795
N 10.62473 -1.18965 1.07813
C 11.41225 -2.10145 0.34157
C 11.06514 -3.45736 0.81427
C 10.07528 -3.24003 1.78454
C 9.80883 -1.84840 1.89350
C 11.60696 -4.67592 0.31660
C 9.18762 -3.87219 2.73381
O 9.02965 -5.02163 3.14318
C 8.13781 -2.79186 3.24095
C 8.06496 -2.88323 4.74057
O 9.04958 -2.58760 5.44037
O 6.78764 -3.09308 5.10770
C 6.66257 -3.17991 6.58711
H 7.14408 3.83633 1.43855
H 12.48135 3.68075 -1.80196
H 12.85729 -2.53499 -0.99089
H 5.98641 -0.51931 3.45410
H 6.00056 2.28327 3.16700
H 5.31956 0.05628 1.29166
H 5.91518 1.53432 0.50748
H 4.49612 1.64729 1.50404
H 7.51052 1.65606 4.90691
H 7.85266 -0.05935 5.33760
H 5.03212 0.32347 5.41199
H 5.69879 1.83156 5.92781
H 7.70135 5.86860 1.42990
H 8.43697 6.76813 0.14580
H 9.18854 6.85480 1.75652
H 10.34278 7.74621 -0.47598
H 12.15293 8.15532 -0.52783
H 11.34828 7.53566 0.94253
H 14.23255 1.90070 -2.31610
H 13.31903 -0.84000 -3.03083
H 13.69759 2.15256 -4.50231
H 11.97624 1.93498 -4.02731
H 12.87725 0.59384 -4.52243
H 15.73110 -0.20776 -2.63841
H 15.21422 -1.68497 -1.84742
H 15.80277 1.00035 -0.83509
H 16.33873 -0.57470 -0.14460
H 14.65960 0.06952 0.13793
H 11.81908 -5.29666 1.17695
H 12.50109 -4.50000 -0.27783
H 10.77049 -5.02632 -0.29397
H 7.18062 -3.01893 2.76999
H 6.17862 -4.13554 6.81128
H 6.13014 -2.26088 6.89636
H 7.57650 -3.26541 7.16319

