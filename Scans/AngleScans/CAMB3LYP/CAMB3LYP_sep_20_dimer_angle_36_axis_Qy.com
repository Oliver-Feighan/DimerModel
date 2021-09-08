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
Mg 10.52029 0.95083 0.80000
C 11.88684 -1.17087 3.35146
C 10.81956 3.53911 3.00134
C 8.90916 2.75398 -1.25435
C 9.28978 -1.92457 -0.73741
N 11.27115 1.04635 2.87830
C 11.82239 0.18818 3.72330
C 12.40617 0.83527 5.00783
C 11.51107 2.16831 5.03369
C 11.20557 2.28657 3.54193
C 10.19740 1.99677 5.77813
C 13.93249 1.17090 4.86740
C 14.63822 1.51732 6.18029
H 16.02560 1.09717 6.35665
N 10.38598 2.92556 0.65506
C 10.56023 3.82774 1.64348
C 10.33318 5.14545 1.03008
C 9.84924 4.93780 -0.36183
C 9.67232 3.49948 -0.36758
C 10.65005 6.43466 1.73076
C 9.57698 5.87949 -1.49554
O 8.93497 5.50736 -2.54060
C 10.00221 7.27403 -1.30598
N 9.10107 0.50804 -0.65330
C 8.59825 1.47165 -1.44485
C 7.80909 0.88120 -2.58895
C 8.08993 -0.69236 -2.53764
C 8.95973 -0.75274 -1.26281
C 6.30134 1.13016 -2.55549
C 8.63235 -1.35310 -3.83451
C 10.03101 -0.85600 -4.12093
N 10.71654 -1.17956 1.10159
C 10.09530 -2.20706 0.35837
C 10.44090 -3.48649 1.01109
C 11.22651 -3.10991 2.11086
C 11.34429 -1.69433 2.15297
C 9.99761 -4.77566 0.60164
C 11.97414 -3.58737 3.25182
O 12.30001 -4.69582 3.67478
C 12.32998 -2.34585 4.17831
C 13.78715 -2.43852 4.53969
O 14.66372 -2.32626 3.66468
O 13.91051 -2.42625 5.87946
C 15.33865 -2.50515 6.28692
H 10.65848 4.37224 3.68604
H 8.42234 3.32196 -2.05052
H 9.02611 -2.87215 -1.21199
H 12.24963 0.26463 5.91538
H 12.08755 3.02463 5.37390
H 10.02883 0.96983 6.06645
H 9.42667 2.32939 5.09520
H 10.15689 2.67776 6.63308
H 14.04182 2.12596 4.33763
H 14.45446 0.37284 4.37620
H 14.03862 1.24069 7.05603
H 14.72829 2.60411 6.24959
H 10.83501 6.27623 2.80038
H 9.74227 7.04690 1.69985
H 11.46372 6.98469 1.25318
H 9.51336 7.68412 -0.42546
H 9.80366 7.77194 -2.25705
H 11.07888 7.28696 -1.10291
H 8.15584 1.26926 -3.53883
H 7.17475 -1.26274 -2.32096
H 5.92170 1.63296 -3.46495
H 6.07126 1.71392 -1.67411
H 5.68964 0.24140 -2.41011
H 8.01888 -1.06471 -4.69162
H 8.64186 -2.43725 -3.79550
H 9.85496 0.09381 -4.63468
H 10.56356 -1.55756 -4.75988
H 10.56766 -0.63343 -3.19134
H 10.85522 -5.43268 0.65944
H 9.58079 -4.75220 -0.40309
H 9.23298 -4.96859 1.35886
H 11.68614 -2.39774 5.05715
H 15.43226 -3.36404 6.95866
H 15.58586 -1.51060 6.70346
H 16.06496 -2.75456 5.52225

