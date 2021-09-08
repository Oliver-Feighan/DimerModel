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
Mg 13.13487 1.18715 1.00700
C 10.89927 -1.38722 2.12675
C 10.44015 3.27590 0.90603
C 14.98479 3.34395 -0.40487
C 15.36834 -1.33428 0.11274
N 10.99763 0.86899 1.47951
C 10.27507 -0.13410 1.95498
C 8.82026 0.24688 2.33960
C 8.64903 1.53093 1.39052
C 10.11421 1.94617 1.27438
C 8.13809 1.18837 0.00094
C 8.68354 0.63033 3.85496
C 7.24495 0.71332 4.37005
H 6.97082 0.28993 5.74047
N 12.89834 3.14723 0.80521
C 11.73505 3.83133 0.81088
C 12.09716 5.24245 0.60649
C 13.55474 5.31523 0.31515
C 13.87210 3.90560 0.20423
C 11.11640 6.36659 0.77332
C 14.50314 6.46384 0.14949
O 15.67838 6.30206 -0.33577
C 13.98269 7.79735 0.48471
N 14.82522 1.02623 -0.19264
C 15.46828 2.12243 -0.63168
C 16.80317 1.76467 -1.24033
C 17.03049 0.21888 -0.89908
C 15.68815 -0.08568 -0.19855
C 16.92417 1.97557 -2.74939
C 18.34814 -0.15357 -0.16576
C 18.33798 0.41554 1.23474
N 13.24218 -0.95623 1.25691
C 14.25289 -1.82077 0.78242
C 13.82897 -3.20042 1.09751
C 12.58039 -3.04299 1.71781
C 12.24028 -1.66406 1.76645
C 14.54415 -4.38619 0.76800
C 11.46793 -3.73053 2.33295
O 11.23541 -4.89623 2.65047
C 10.27232 -2.69852 2.51148
C 9.74368 -2.83799 3.91276
O 10.45214 -2.53428 4.88870
O 8.42442 -3.09755 3.86222
C 7.85222 -3.23345 5.22822
H 9.63063 3.94870 0.62157
H 15.70938 4.05414 -0.80995
H 16.05413 -2.16933 -0.04557
H 8.07488 -0.49972 2.09311
H 8.06972 2.30949 1.88003
H 8.08658 0.12144 -0.15770
H 8.83817 1.64018 -0.68957
H 7.17755 1.67955 -0.17921
H 8.99172 1.67521 3.98912
H 9.24957 -0.04146 4.47053
H 6.53133 0.25336 3.67569
H 6.94800 1.76442 4.40230
H 10.08517 5.99621 0.82609
H 11.14625 6.95653 -0.14914
H 11.36003 7.01649 1.61649
H 13.11204 8.01045 -0.13093
H 14.83275 8.47555 0.38823
H 13.63785 7.78704 1.52471
H 17.60053 2.34450 -0.79187
H 17.05780 -0.39939 -1.80838
H 17.75730 2.64696 -3.03086
H 15.98303 2.36293 -3.11646
H 17.04349 1.06547 -3.33479
H 19.20491 0.29323 -0.67612
H 18.52602 -1.22250 -0.11340
H 18.67015 1.44710 1.08434
H 19.02665 -0.13144 1.87551
H 17.31895 0.45271 1.63707
H 14.50393 -5.02631 1.63929
H 15.57069 -4.16515 0.48316
H 13.95105 -4.74286 -0.07833
H 9.51667 -2.93992 1.76294
H 7.35965 -4.20971 5.27261
H 7.21558 -2.34086 5.37559
H 8.54627 -3.30886 6.05711

