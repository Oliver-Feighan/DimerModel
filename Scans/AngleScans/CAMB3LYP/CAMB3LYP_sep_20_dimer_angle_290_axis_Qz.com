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
Mg 10.51751 0.95443 0.79724
C 9.91578 2.77077 3.83334
C 9.51187 3.74108 -0.89328
C 10.60540 -0.75550 -1.87406
C 10.36872 -1.97287 2.68257
N 9.82723 2.93774 1.49169
C 9.69822 3.53402 2.66734
C 9.39270 5.05384 2.58760
C 8.73568 5.09051 1.12273
C 9.42302 3.86903 0.51591
C 7.23666 4.84106 1.12907
C 10.69304 5.92931 2.65599
C 10.45183 7.42155 2.89435
H 11.41249 8.16048 3.70895
N 10.59283 1.53679 -1.09847
C 10.16058 2.70966 -1.60721
C 10.42220 2.64177 -3.05328
C 10.89273 1.26845 -3.38087
C 10.69098 0.60957 -2.10601
C 10.29755 3.83199 -3.95951
C 11.42297 0.64728 -4.63753
O 11.54752 -0.62245 -4.76039
C 11.70403 1.56255 -5.75323
N 10.28097 -1.07522 0.41393
C 10.41983 -1.56447 -0.83077
C 10.47999 -3.07333 -0.82023
C 10.59607 -3.48792 0.71998
C 10.48838 -2.08974 1.36721
C 9.28732 -3.79108 -1.45170
C 11.79259 -4.39567 1.11621
C 13.09037 -3.63614 0.96048
N 10.38684 0.46154 2.89678
C 10.31094 -0.82806 3.46708
C 10.11985 -0.65050 4.92117
C 10.07441 0.74106 5.09432
C 10.20621 1.38499 3.83437
C 9.96533 -1.69681 5.87379
C 9.93427 1.81151 6.05516
O 9.87755 1.87521 7.28254
C 9.69850 3.17106 5.26623
C 10.61474 4.21155 5.84967
O 11.84860 4.09891 5.74236
O 9.88437 5.27041 6.24433
C 10.74999 6.33396 6.81997
H 9.00086 4.47732 -1.51433
H 10.74142 -1.39226 -2.75128
H 10.41120 -2.83508 3.35150
H 8.68275 5.41365 3.32264
H 9.01754 5.99257 0.58588
H 6.87877 4.56223 2.10903
H 7.06898 4.03866 0.42261
H 6.70809 5.71076 0.72836
H 11.15028 5.96980 1.65898
H 11.36938 5.54971 3.39700
H 9.42734 7.62522 3.22884
H 10.52991 7.94291 1.93721
H 9.77456 4.66039 -3.46590
H 9.63237 3.54110 -4.77974
H 11.25671 4.14600 -4.37675
H 10.79523 2.09691 -6.01968
H 12.15879 0.94509 -6.53028
H 12.42428 2.31578 -5.41488
H 11.35697 -3.42909 -1.34696
H 9.71749 -4.05816 1.05577
H 9.57445 -4.47677 -2.27113
H 8.58913 -3.04499 -1.80725
H 8.67630 -4.36809 -0.75967
H 11.85392 -5.25583 0.44501
H 11.72588 -4.78096 2.12819
H 13.31347 -3.74464 -0.10514
H 13.87076 -4.08122 1.57443
H 12.94755 -2.56780 1.16025
H 10.57245 -1.43136 6.72903
H 10.26000 -2.65863 5.45925
H 8.88838 -1.63974 6.05322
H 8.64774 3.43932 5.38252
H 10.36877 6.54799 7.82324
H 10.73261 7.16248 6.08716
H 11.78103 6.07904 7.03542

