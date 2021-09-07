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
Mg 15.77329 1.41641 1.20503
C 15.93012 -2.11663 0.59501
C 14.71719 0.67327 4.36209
C 15.19370 4.50663 1.63167
C 15.72577 1.87030 -2.25003
N 15.41003 -0.51041 2.22859
C 15.56277 -1.78859 1.91670
C 15.38557 -2.76711 3.10862
C 14.46612 -1.84152 4.04497
C 14.91180 -0.46891 3.54498
C 12.97982 -1.99189 3.76629
C 16.74829 -3.11419 3.80470
C 16.69240 -4.29709 4.77392
H 17.85380 -5.17760 4.86520
N 15.53614 2.45980 2.87646
C 15.12170 1.99524 4.07383
C 15.09428 3.16138 4.97023
C 15.36789 4.38074 4.16213
C 15.35567 3.81972 2.82602
C 14.90589 3.04234 6.45474
C 15.58953 5.81662 4.52998
O 15.57354 6.74757 3.64912
C 15.73282 6.11163 5.96314
N 15.28615 2.95387 -0.10693
C 15.13944 4.21644 0.33160
C 15.06457 5.18905 -0.82114
C 15.43577 4.33991 -2.12450
C 15.58024 2.94228 -1.48331
C 13.71340 5.87323 -1.02757
C 16.59946 4.87342 -3.00420
C 17.90633 4.76592 -2.25192
N 16.00344 0.15362 -0.53291
C 15.92253 0.54391 -1.88761
C 16.03028 -0.68415 -2.70168
C 16.14419 -1.71431 -1.75601
C 16.08956 -1.16899 -0.44485
C 15.96766 -0.74868 -4.12225
C 16.28754 -3.14017 -1.56870
O 16.47295 -4.10013 -2.31558
C 16.02825 -3.47726 -0.03719
C 17.12907 -4.38939 0.43049
O 18.30430 -3.98611 0.48356
O 16.58606 -5.50956 0.94111
C 17.63799 -6.44105 1.42854
H 14.16500 0.55585 5.29492
H 15.10134 5.59155 1.72050
H 15.81789 1.93758 -3.33615
H 14.86486 -3.68626 2.86821
H 14.71987 -1.96588 5.09449
H 12.79269 -2.61056 2.90121
H 12.60846 -0.98889 3.60160
H 12.46468 -2.36144 4.65757
H 17.00783 -2.30972 4.50485
H 17.51713 -3.28057 3.07543
H 15.77351 -4.88405 4.65536
H 16.63225 -3.90799 5.79323
H 14.56380 2.03918 6.73764
H 14.07249 3.69979 6.72504
H 15.79163 3.34501 7.01734
H 14.84161 5.78042 6.49088
H 15.97290 7.17512 6.02001
H 16.56931 5.52506 6.35914
H 15.78879 5.98565 -0.70267
H 14.59013 4.27832 -2.82517
H 13.77538 6.97764 -1.00478
H 13.03000 5.51501 -0.26916
H 13.19689 5.60716 -1.94824
H 16.45348 5.93412 -3.22272
H 16.69347 4.35878 -3.95462
H 17.91048 5.66896 -1.63418
H 18.74663 4.75672 -2.94319
H 17.90108 3.90360 -1.57514
H 16.74873 -1.42754 -4.43781
H 16.08974 0.23506 -4.57070
H 14.95460 -1.13762 -4.25535
H 15.04768 -3.94967 0.03281
H 17.47684 -7.39794 0.92250
H 17.55514 -6.43436 2.53164
H 18.66403 -6.22650 1.15346

