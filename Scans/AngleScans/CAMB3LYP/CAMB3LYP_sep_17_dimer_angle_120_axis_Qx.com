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
Mg 8.92656 0.80880 0.68447
C 10.85426 2.30951 3.31333
C 7.33089 -1.00453 3.09285
C 7.53195 -0.85302 -1.63087
C 11.39893 1.85596 -1.53736
N 9.15139 0.74907 2.88405
C 9.87319 1.40869 3.77750
C 9.46900 1.14431 5.25273
C 8.78511 -0.29587 5.06028
C 8.36084 -0.17160 3.59836
C 9.76796 -1.44911 5.17600
C 8.44098 2.20126 5.78942
C 8.23872 2.18903 7.30622
H 8.00112 3.45890 7.98699
N 7.30515 -0.33496 0.72210
C 6.78591 -0.97565 1.79036
C 5.60561 -1.69829 1.29151
C 5.56546 -1.56482 -0.19015
C 6.83857 -0.91575 -0.43094
C 4.60923 -2.35570 2.20170
C 4.55495 -1.96771 -1.22105
O 4.83301 -1.96837 -2.47225
C 3.26096 -2.45049 -0.71707
N 9.49462 0.35314 -1.26284
C 8.67972 -0.34145 -2.07611
C 9.17587 -0.30992 -3.50203
C 10.36742 0.75662 -3.51884
C 10.41483 1.10679 -2.01533
C 9.67685 -1.64117 -4.06164
C 10.25842 1.91750 -4.54506
C 9.10974 2.82778 -4.17452
N 10.69755 2.04134 0.79637
C 11.60422 2.32856 -0.24744
C 12.71028 3.11216 0.33994
C 12.37930 3.20756 1.70003
C 11.16379 2.51395 1.94690
C 13.85632 3.58638 -0.35851
C 12.78723 3.72826 2.98500
O 13.68859 4.46895 3.37573
C 11.88314 3.06260 4.11013
C 11.39960 4.15504 5.02407
O 10.62708 5.03458 4.60426
O 11.75217 3.85991 6.28857
C 11.28507 4.90873 7.23388
H 6.92879 -1.78375 3.74093
H 7.03639 -1.31734 -2.48662
H 12.15818 2.30334 -2.18239
H 10.29908 1.07220 5.94522
H 7.91486 -0.40875 5.70160
H 10.78809 -1.10303 5.25142
H 9.63437 -2.04054 4.27967
H 9.48966 -2.09930 6.01035
H 7.43220 1.91736 5.46306
H 8.70778 3.18938 5.46840
H 9.01976 1.61694 7.82161
H 7.31727 1.64562 7.52907
H 4.99027 -2.43134 3.22761
H 4.51133 -3.39536 1.87097
H 3.62717 -1.87934 2.16541
H 3.42265 -3.31241 -0.07410
H 2.64173 -2.60528 -1.60277
H 2.81612 -1.66943 -0.09042
H 8.39604 0.02376 -4.17555
H 11.33153 0.28662 -3.76259
H 9.15031 -1.95300 -4.98338
H 9.57512 -2.39485 -3.29213
H 10.74243 -1.67795 -4.28161
H 10.03425 1.52401 -5.53958
H 11.16558 2.50661 -4.62877
H 8.24226 2.31775 -4.60398
H 9.24531 3.81394 -4.61398
H 8.96828 2.86059 -3.08798
H 14.02297 4.60334 -0.02909
H 13.71475 3.53486 -1.43604
H 14.60710 2.87285 -0.00855
H 12.49918 2.33193 4.63566
H 12.16473 5.25876 7.78283
H 10.47629 4.43941 7.82495
H 10.91661 5.83617 6.81133

