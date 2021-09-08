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
Mg 15.76230 1.42735 1.20343
C 18.00440 1.41994 4.00557
C 13.23493 0.82955 3.41469
C 13.83270 0.98063 -1.27538
C 18.53383 0.88404 -0.83841
N 15.74716 1.17797 3.40112
C 16.65066 1.25688 4.36643
C 16.06129 1.23581 5.80231
C 14.67915 0.48269 5.48391
C 14.51710 0.87875 4.01777
C 14.78713 -1.03117 5.56058
C 15.81250 2.67825 6.36754
C 15.53084 2.74032 7.87039
H 16.03394 3.88406 8.62617
N 13.78086 1.45908 1.09626
C 12.90786 1.21482 2.09601
C 11.56860 1.34308 1.50095
C 11.72324 1.52068 0.03151
C 13.14852 1.30140 -0.11196
C 10.31251 1.37158 2.32245
C 10.74752 1.82349 -1.06506
O 11.06326 1.69804 -2.30094
C 9.38360 2.18109 -0.64871
N 16.09623 0.78540 -0.74560
C 15.09124 0.73066 -1.63719
C 15.61461 0.50804 -3.03609
C 17.20249 0.66679 -2.93105
C 17.33574 0.87456 -1.40643
C 15.27908 -0.84432 -3.66452
C 17.87071 1.70000 -3.87902
C 17.45077 3.10039 -3.49440
N 17.90547 1.37441 1.47185
C 18.88190 1.10420 0.48815
C 20.19054 1.06616 1.17258
C 19.87966 1.29588 2.52131
C 18.47428 1.44506 2.67015
C 21.44436 0.79495 0.55575
C 20.41976 1.43552 3.85454
O 21.55254 1.48967 4.33163
C 19.21751 1.39625 4.89350
C 19.40105 2.53453 5.85950
O 19.32454 3.71297 5.46959
O 19.41955 2.04901 7.11421
C 19.58847 3.14142 8.10916
H 12.40602 0.41753 3.99091
H 13.22347 0.92477 -2.18048
H 19.45599 0.81758 -1.41967
H 16.63775 0.66658 6.52165
H 13.86389 0.88470 6.07982
H 15.80745 -1.35570 5.70121
H 14.39832 -1.40178 4.62115
H 14.12009 -1.41786 6.33637
H 14.85624 3.05409 5.98148
H 16.63078 3.32844 6.12609
H 15.78753 1.80200 8.37697
H 14.45251 2.83809 8.01760
H 10.50021 1.05386 3.35552
H 9.64810 0.60053 1.91761
H 9.80238 2.33621 2.27867
H 8.96135 1.37005 -0.06017
H 8.85834 2.44890 -1.56755
H 9.43684 3.05420 0.01118
H 15.23070 1.25816 -3.71656
H 17.72244 -0.27339 -3.16672
H 14.73913 -0.75678 -4.62622
H 14.69854 -1.41662 -2.95320
H 16.13321 -1.49542 -3.84266
H 17.53212 1.54562 -4.90649
H 18.95410 1.64319 -3.88499
H 16.48343 3.21363 -3.99292
H 18.17073 3.83008 -3.85942
H 17.27615 3.17638 -2.41487
H 22.15102 1.50719 0.96055
H 21.37904 0.87045 -0.52768
H 21.60582 -0.23451 0.88650
H 19.24736 0.42692 5.39258
H 20.46303 2.88813 8.71634
H 18.61745 3.22107 8.63311
H 19.86617 4.12037 7.73634

