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
Mg 15.77557 1.42402 1.19983
C 16.99423 -0.70195 3.82169
C 15.91752 4.00247 3.42830
C 14.28385 3.22712 -0.94284
C 14.66458 -1.45142 -0.42597
N 16.39370 1.51353 3.32166
C 16.89687 0.65462 4.19541
C 17.39419 1.29880 5.51726
C 16.48975 2.62556 5.49335
C 16.27773 2.74971 3.98593
C 15.13317 2.44111 6.15290
C 18.92388 1.64555 5.47467
C 19.54327 1.98974 6.83101
H 20.91979 1.57810 7.09216
N 15.63656 3.39853 1.05643
C 15.74195 4.29659 2.05827
C 15.54449 5.61598 1.43828
C 15.15043 5.41249 0.01771
C 14.98447 3.97305 -0.00621
C 15.80753 6.90355 2.16382
C 14.94322 6.35835 -1.12623
O 14.37078 5.98745 -2.21139
C 15.34575 7.75472 -0.90348
N 14.45363 0.97936 -0.34195
C 13.99467 1.94375 -1.15879
C 13.28315 1.35406 -2.35312
C 13.57143 -0.21780 -2.29198
C 14.35984 -0.27907 -0.96532
C 11.77455 1.59258 -2.41326
C 14.19893 -0.86789 -3.55542
C 15.60924 -0.35976 -3.75093
N 15.96767 -0.70657 1.50270
C 15.40167 -1.73429 0.71687
C 15.71470 -3.01480 1.38374
C 16.42697 -2.63877 2.53254
C 16.53178 -1.22267 2.58892
C 15.30720 -4.30475 0.94091
C 17.10484 -3.11723 3.71587
O 17.41139 -4.22568 4.15303
C 17.39292 -1.87829 4.66898
C 18.82514 -1.96297 5.12076
O 19.75412 -1.84007 4.30314
O 18.86401 -1.95703 6.46569
C 20.26426 -2.02839 6.96170
H 15.70781 4.83081 4.10561
H 13.84397 3.79603 -1.76523
H 14.43799 -2.39823 -0.92083
H 17.18503 0.72226 6.41036
H 17.03761 3.48395 5.87332
H 14.95415 1.41152 6.42502
H 14.40451 2.77212 5.42453
H 15.03419 3.11723 7.00695
H 19.05947 2.60415 4.95752
H 19.48135 0.85370 5.01335
H 18.89183 1.70436 7.66596
H 19.62106 3.07672 6.91117
H 15.92606 6.74067 3.24216
H 14.89914 7.50976 2.07892
H 16.64565 7.46164 1.74102
H 14.79965 8.15676 -0.05342
H 15.20379 8.25634 -1.86269
H 16.40742 7.77388 -0.63308
H 13.68611 1.74955 -3.27741
H 12.64853 -0.79554 -2.13605
H 11.44921 2.09765 -3.34230
H 11.48541 2.17004 -1.54523
H 11.16127 0.69892 -2.31098
H 13.63847 -0.57911 -4.44797
H 14.21369 -1.95214 -3.52122
H 15.45905 0.59157 -4.27006
H 16.18586 -1.05426 -4.35858
H 16.08484 -0.13852 -2.78837
H 16.16415 -4.95622 1.04927
H 14.95415 -4.27876 -0.08790
H 14.49791 -4.50692 1.64762
H 16.69554 -1.93925 5.50535
H 20.32161 -2.89020 7.63378
H 20.47772 -1.03442 7.39783
H 21.03893 -2.26876 6.24298

