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
Mg 10.51480 0.94027 0.80624
C 7.82450 1.05543 3.17860
C 12.49100 2.34723 3.20410
C 12.77858 1.26830 -1.39266
C 8.10985 0.62676 -1.69648
N 10.09701 1.57310 2.88404
C 9.05406 1.51117 3.69814
C 9.36511 1.88731 5.17168
C 10.66219 2.79959 4.91832
C 11.14077 2.18461 3.60477
C 10.33036 4.26136 4.66841
C 9.70787 0.63258 6.04927
C 9.71897 0.88682 7.55824
H 9.25338 -0.17473 8.44639
N 12.46739 1.23299 1.00697
C 13.10203 1.79728 2.05582
C 14.53133 1.79923 1.70760
C 14.67152 1.34061 0.29879
C 13.27925 1.28142 -0.09894
C 15.60945 2.13392 2.69707
C 15.86219 1.02194 -0.55382
O 15.76071 0.87041 -1.82251
C 17.16482 0.98475 0.12681
N 10.45198 1.15999 -1.26019
C 11.58379 1.22906 -1.98279
C 11.29603 1.10699 -3.46018
C 9.75371 0.69734 -3.56506
C 9.37748 0.74256 -2.06780
C 11.54615 2.36424 -4.29274
C 9.41997 -0.58687 -4.37259
C 9.95608 -1.80321 -3.65254
N 8.37263 0.66867 0.73287
C 7.56208 0.58784 -0.42062
C 6.15911 0.51957 0.03703
C 6.25010 0.58917 1.43537
C 7.61217 0.71379 1.82105
C 5.01205 0.45662 -0.80346
C 5.50287 0.59625 2.67242
O 4.32155 0.43238 2.97467
C 6.47973 1.02575 3.85041
C 6.28642 0.06189 4.98888
O 6.59710 -1.13569 4.86280
O 5.97322 0.75098 6.10128
C 5.78245 -0.16080 7.26063
H 13.13573 2.99158 3.80245
H 13.52745 1.26768 -2.18806
H 7.30799 0.42904 -2.41104
H 8.59385 2.46644 5.66519
H 11.40452 2.65631 5.69921
H 9.26612 4.42436 4.58504
H 10.82663 4.52155 3.74263
H 10.78500 4.88848 5.44065
H 10.76137 0.36540 5.89600
H 9.04643 -0.18153 5.82485
H 9.24557 1.84122 7.81894
H 10.75575 1.00446 7.88271
H 15.19523 2.59199 3.60371
H 16.22100 2.92427 2.24850
H 16.24957 1.27965 2.92759
H 17.35537 1.94763 0.59479
H 17.87994 0.65356 -0.62870
H 17.11642 0.24621 0.93480
H 11.89827 0.32619 -3.90825
H 9.15801 1.47682 -4.06245
H 12.25797 2.20269 -5.12417
H 11.90215 3.14486 -3.63363
H 10.65541 2.81603 -4.72594
H 9.91366 -0.56043 -5.34713
H 8.35848 -0.72010 -4.55284
H 11.00377 -1.83392 -3.96636
H 9.42387 -2.69921 -3.96557
H 9.94149 -1.65581 -2.56639
H 4.35112 -0.28404 -0.37325
H 5.28219 0.20320 -1.82648
H 4.65032 1.48471 -0.71709
H 6.22413 2.04807 4.13166
H 4.78573 0.04085 7.66505
H 6.64410 0.02246 7.92977
H 5.71720 -1.22326 7.05725

