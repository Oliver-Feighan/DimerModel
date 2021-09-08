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
Mg 13.12952 1.18517 0.99305
C 12.85966 -1.60172 -1.25184
C 13.17276 3.20383 -1.75608
C 13.81845 3.61640 2.91187
C 14.20089 -1.06161 3.43228
N 13.02287 0.74048 -1.17101
C 12.85968 -0.35086 -1.90360
C 12.60635 -0.08156 -3.41121
C 13.29343 1.36615 -3.51565
C 13.12544 1.82149 -2.06751
C 14.77940 1.30594 -3.82808
C 11.07748 -0.01818 -3.75841
C 10.75499 -0.06691 -5.25353
H 9.54680 -0.76598 -5.68254
N 12.96487 3.14077 0.69795
C 12.97930 3.78918 -0.48554
C 12.84634 5.21831 -0.16294
C 12.91795 5.37413 1.31534
C 13.26497 4.02255 1.70637
C 12.59168 6.26762 -1.20578
C 12.72090 6.54849 2.22542
O 13.06572 6.50949 3.45928
C 12.19883 7.77625 1.60772
N 14.09796 1.30181 2.82903
C 14.21787 2.47547 3.47406
C 14.70014 2.27260 4.89063
C 14.64775 0.69384 5.14036
C 14.20991 0.22956 3.73390
C 16.10605 2.78938 5.19515
C 13.83298 0.20216 6.36799
C 12.36134 0.47291 6.15272
N 13.30158 -0.96190 1.16206
C 13.79881 -1.69376 2.26262
C 13.82637 -3.11577 1.86317
C 13.36127 -3.10987 0.53948
C 13.08390 -1.77601 0.13527
C 14.28706 -4.19674 2.66653
C 13.05195 -3.92761 -0.61131
O 13.00852 -5.13741 -0.83092
C 12.83335 -2.97455 -1.86442
C 11.57269 -3.40640 -2.56211
O 10.46928 -3.29633 -1.99911
O 11.85160 -3.67282 -3.85114
C 10.63035 -4.09472 -4.58760
H 13.42330 3.90653 -2.55128
H 13.97168 4.40672 3.65038
H 14.42919 -1.83479 4.16900
H 13.09552 -0.77478 -4.08495
H 12.74066 2.01820 -4.18695
H 15.15387 0.29320 -3.81151
H 15.26464 1.90125 -3.06572
H 14.98543 1.80697 -4.77827
H 10.69656 0.98240 -3.51656
H 10.53825 -0.79171 -3.24704
H 11.61712 -0.38753 -5.85092
H 10.55333 0.95017 -5.59836
H 12.75367 5.87659 -2.21783
H 13.36475 7.03399 -1.08372
H 11.61048 6.73567 -1.10153
H 12.87102 8.09561 0.81485
H 12.04207 8.47438 2.43229
H 11.23815 7.54777 1.13285
H 14.04134 2.76222 5.59729
H 15.65025 0.27493 5.31164
H 16.13245 3.51834 6.02702
H 16.51163 3.22731 4.29279
H 16.84023 2.02025 5.42887
H 14.12329 0.75847 7.26264
H 13.97216 -0.85170 6.58477
H 12.25600 1.51900 6.45551
H 11.75614 -0.18150 6.77683
H 12.10088 0.40957 5.08979
H 13.57509 -5.00194 2.54300
H 14.38081 -3.90455 3.71037
H 15.25975 -4.38314 2.20349
H 13.70997 -3.07090 -2.50616
H 10.84912 -5.06778 -5.03812
H 10.39451 -3.26215 -5.27677
H 9.74747 -4.32593 -4.00321

