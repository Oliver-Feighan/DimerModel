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
Mg 7.89077 0.71193 0.59182
C 7.39884 4.23975 0.15417
C 7.29251 0.36552 -2.74834
C 7.87211 -2.43161 1.02222
C 7.31957 1.19531 3.99561
N 7.41901 2.18673 -0.98763
C 7.31878 3.50607 -1.04780
C 7.20188 4.07729 -2.48641
C 6.60952 2.77455 -3.21501
C 7.16690 1.69936 -2.28445
C 5.09218 2.70057 -3.17556
C 8.59478 4.48538 -3.08250
C 8.52642 5.33255 -4.35505
H 9.53788 6.36603 -4.55837
N 8.10683 -0.76472 -0.71637
C 7.83850 -0.72638 -2.03842
C 8.14412 -2.07133 -2.55015
C 8.45765 -2.95335 -1.39323
C 8.13025 -2.07005 -0.29201
C 8.19794 -2.39150 -4.01576
C 8.95593 -4.36337 -1.29439
O 8.93105 -5.01032 -0.18819
C 9.38694 -4.99588 -2.54964
N 7.42587 -0.48454 2.22726
C 7.55137 -1.82174 2.16349
C 7.42878 -2.44704 3.53257
C 7.43615 -1.22567 4.56494
C 7.47356 -0.04926 3.56490
C 6.18284 -3.30157 3.76447
C 8.49882 -1.26749 5.69705
C 9.88361 -1.11141 5.11107
N 7.62294 2.41203 1.89796
C 7.37080 2.38900 3.28720
C 7.14807 3.78613 3.71247
C 7.26423 4.52857 2.52766
C 7.51988 3.65368 1.43738
C 6.83293 4.21351 5.03313
C 7.21818 5.85239 1.94966
O 7.12352 6.99374 2.39922
C 7.17777 5.71165 0.36696
C 8.19020 6.66280 -0.20994
O 9.40463 6.49749 0.00083
O 7.57848 7.47791 -1.08845
C 8.54377 8.43252 -1.69570
H 6.89764 0.11985 -3.73455
H 7.96321 -3.49730 1.24464
H 7.23269 1.43933 5.05654
H 6.51340 4.90723 -2.59108
H 7.01799 2.66169 -4.21590
H 4.66677 3.47857 -2.55920
H 4.85564 1.72531 -2.77065
H 4.68771 2.71119 -4.19174
H 9.09066 3.59013 -3.47920
H 9.19214 4.98695 -2.34615
H 7.52212 5.73873 -4.52637
H 8.70263 4.68183 -5.21501
H 7.76046 -1.58689 -4.61960
H 7.53380 -3.24680 -4.18081
H 9.20293 -2.65421 -4.35289
H 8.55982 -4.99797 -3.25556
H 9.79280 -5.96791 -2.26275
H 10.17952 -4.38449 -2.99532
H 8.27563 -3.08907 3.74171
H 6.48327 -1.14586 5.10851
H 6.41577 -4.33791 4.07401
H 5.59350 -3.29908 2.85716
H 5.48119 -2.90586 4.49679
H 8.48117 -2.23886 6.19732
H 8.34822 -0.51004 6.45900
H 10.13183 -2.12759 4.79036
H 10.58107 -0.75437 5.86614
H 9.86320 -0.47897 4.21594
H 7.43509 5.09086 5.22847
H 7.02486 3.42840 5.76151
H 5.76474 4.41828 4.92218
H 6.16339 5.94741 0.04289
H 8.15374 9.43910 -1.51579
H 8.65389 8.11596 -2.74994
H 9.52775 8.50368 -1.24718

