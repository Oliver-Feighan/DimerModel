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
Mg 7.89049 0.70445 0.60437
C 8.10022 -2.34727 -1.27236
C 6.97339 -1.18599 3.29133
C 7.26095 3.39787 2.15890
C 7.69982 2.39808 -2.43552
N 7.60887 -1.47526 0.85302
C 7.77713 -2.54239 0.08668
C 7.66732 -3.89684 0.83698
C 6.76484 -3.41382 2.07435
C 7.16104 -1.93911 2.10486
C 5.27251 -3.50190 1.80136
C 9.06310 -4.42898 1.31724
C 9.07031 -5.88702 1.78179
H 10.25318 -6.69866 1.50877
N 7.69546 1.04940 2.54899
C 7.33772 0.16222 3.50093
C 7.31964 0.91411 4.76525
C 7.53537 2.35412 4.45775
C 7.48375 2.32508 3.00979
C 7.19132 0.24982 6.10527
C 7.73986 3.55999 5.32392
O 7.66971 4.74895 4.85032
C 7.93192 3.31067 6.76004
N 7.31987 2.59953 -0.03234
C 7.16285 3.60568 0.84569
C 7.02248 4.93153 0.13661
C 7.36142 4.63595 -1.39810
C 7.56079 3.10624 -1.32316
C 5.64983 5.59671 0.23540
C 8.47845 5.49586 -2.05039
C 9.81540 5.16365 -1.42783
N 8.08071 0.17996 -1.48292
C 7.93913 1.03911 -2.59463
C 8.04198 0.20211 -3.80759
C 8.21461 -1.09953 -3.31320
C 8.19882 -1.07816 -1.89236
C 7.92599 0.66367 -5.14910
C 8.39596 -2.48818 -3.67021
O 8.57317 -3.09823 -4.72386
C 8.20334 -3.37479 -2.36514
C 9.34091 -4.35664 -2.29857
O 10.50829 -3.96116 -2.13297
O 8.84242 -5.60412 -2.22319
C 9.93233 -6.61291 -2.14406
H 6.46029 -1.65793 4.12972
H 7.14859 4.36978 2.64496
H 7.74848 2.86412 -3.42199
H 7.15777 -4.67997 0.28849
H 7.06155 -3.90750 2.99609
H 5.06558 -3.76415 0.77445
H 4.87341 -2.52224 2.02945
H 4.80028 -4.19147 2.50693
H 9.33196 -3.93084 2.25762
H 9.80664 -4.28826 0.55699
H 8.16047 -6.42018 1.48034
H 9.04113 -5.90339 2.87398
H 6.88223 -0.79812 6.00679
H 6.35491 0.73224 6.62253
H 8.09134 0.35406 6.71499
H 7.06914 2.77780 7.15261
H 8.15091 4.28588 7.19922
H 8.79560 2.64857 6.88774
H 7.73329 5.65276 0.52088
H 6.49090 4.80793 -2.04807
H 5.68870 6.61637 0.66302
H 5.00414 4.96086 0.82643
H 5.10407 5.67115 -0.70364
H 8.30118 6.55672 -1.85725
H 8.54683 5.37138 -3.12594
H 9.82382 5.77494 -0.52045
H 10.62839 5.43882 -2.09668
H 9.85497 4.11294 -1.11781
H 8.70883 0.17624 -5.71484
H 8.00932 1.74702 -5.20542
H 6.91722 0.31658 -5.38824
H 7.23668 -3.87317 -2.44741
H 9.77255 -7.32081 -2.96322
H 9.89199 -7.01618 -1.11474
H 10.94209 -6.27691 -2.34884

