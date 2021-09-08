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
Mg 31.53887 2.84411 2.40198
C 31.29150 0.80381 5.34393
C 30.28317 5.41407 4.26024
C 31.32039 4.53907 -0.27137
C 31.67611 -0.09929 0.54100
N 30.90786 2.96946 4.51820
C 30.91702 2.14697 5.55635
C 30.57856 2.82033 6.91341
C 29.74289 4.07411 6.35793
C 30.36922 4.18962 4.96976
C 28.26082 3.78293 6.19051
C 31.86242 3.28332 7.68759
C 31.62803 3.66673 9.15040
H 32.67058 3.36025 10.12591
N 31.36665 4.80292 2.13373
C 30.85460 5.70111 3.00111
C 30.93431 7.00478 2.32416
C 31.38304 6.77923 0.92327
C 31.35344 5.33146 0.86686
C 30.67737 8.30486 3.02928
C 31.76217 7.71249 -0.18646
O 31.89158 7.30461 -1.39464
C 31.88608 9.13436 0.16625
N 31.29019 2.30774 0.40879
C 31.28037 3.23848 -0.56154
C 31.36659 2.60240 -1.92842
C 31.68140 1.05690 -1.66399
C 31.63218 1.04913 -0.12033
C 30.11645 2.72510 -2.79932
C 32.93310 0.46339 -2.36661
C 34.18832 1.07896 -1.79138
N 31.68020 0.73289 2.83884
C 31.69782 -0.33096 1.91040
C 31.68724 -1.58727 2.68757
C 31.64058 -1.16391 4.02438
C 31.60291 0.25574 4.07610
C 31.67066 -2.90227 2.14292
C 31.60389 -1.59534 5.40323
O 31.70014 -2.67937 5.97712
C 31.24856 -0.34302 6.31529
C 32.21309 -0.32026 7.46928
O 33.42553 -0.12234 7.27583
O 31.52024 -0.31846 8.62268
C 32.43319 -0.28509 9.79624
H 29.68623 6.22036 4.68737
H 31.35620 5.07697 -1.22156
H 31.81637 -1.05804 0.03721
H 29.96026 2.22478 7.57431
H 29.94024 4.96989 6.94103
H 28.03308 2.74124 6.36131
H 28.01950 4.06106 5.17294
H 27.67180 4.44498 6.83166
H 32.19057 4.25220 7.28967
H 32.63242 2.53888 7.62807
H 30.65019 3.32795 9.51371
H 31.57700 4.75605 9.21972
H 30.20648 8.14584 4.00723
H 29.91974 8.84023 2.44685
H 31.57139 8.92766 3.10347
H 30.94228 9.48896 0.57363
H 32.25029 9.62871 -0.73644
H 32.63032 9.23328 0.96434
H 32.17659 3.03242 -2.50469
H 30.86139 0.40942 -2.00781
H 30.31020 3.20869 -3.77537
H 29.36688 3.27416 -2.24505
H 29.60656 1.78661 -3.01015
H 32.92433 0.71315 -3.43042
H 32.99999 -0.61672 -2.28994
H 34.27810 2.01873 -2.34442
H 35.04672 0.43383 -1.96742
H 34.05377 1.33086 -0.73304
H 32.37550 -3.48682 2.71911
H 31.92717 -2.89794 1.08559
H 30.62382 -3.16949 2.31025
H 30.21628 -0.46158 6.64681
H 32.17698 -1.14351 10.42484
H 32.31040 0.71999 10.24177
H 33.48720 -0.46169 9.61646

