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
Mg 10.51891 0.94286 0.79246
C 10.31699 3.21664 -1.97669
C 10.03333 -1.58752 -1.44271
C 10.24191 -1.04289 3.25163
C 9.84734 3.64859 2.88381
N 10.22665 0.94105 -1.39999
C 10.23032 1.86024 -2.35363
C 10.21660 1.28950 -3.79711
C 9.55657 -0.14059 -3.48399
C 9.99134 -0.29699 -2.02825
C 8.03768 -0.12542 -3.52984
C 11.65992 1.13778 -4.39361
C 11.70888 0.88015 -5.90125
H 12.80402 1.46283 -6.67173
N 10.67470 -1.03405 0.87068
C 10.46452 -1.90735 -0.13652
C 10.68695 -3.24378 0.43675
C 10.88430 -3.09766 1.90454
C 10.58065 -1.69066 2.07257
C 10.77614 -4.48486 -0.40301
C 11.26865 -4.06705 2.98100
O 11.14896 -3.77584 4.22346
C 11.70105 -5.40069 2.53838
N 9.89708 1.21097 2.75811
C 9.92229 0.19296 3.63615
C 9.69615 0.68324 5.04642
C 9.75471 2.27913 4.96103
C 9.92308 2.44488 3.43474
C 8.37996 0.25656 5.69601
C 10.76373 2.99751 5.89830
C 12.17927 2.67000 5.48072
N 10.32876 3.08199 0.55577
C 10.01885 4.02690 1.55835
C 9.88653 5.33952 0.89355
C 10.10770 5.06108 -0.46367
C 10.34007 3.66968 -0.63543
C 9.55116 6.56601 1.53336
C 10.18691 5.62616 -1.79149
O 10.16160 6.76626 -2.25331
C 10.20077 4.43747 -2.84657
C 11.30587 4.70363 -3.83167
O 12.49441 4.69503 -3.46575
O 10.79492 4.70847 -5.07637
C 11.85457 4.95759 -6.08974
H 9.66157 -2.43272 -2.02264
H 10.24192 -1.66618 4.14884
H 9.73598 4.55720 3.47937
H 9.59857 1.83902 -4.49701
H 9.99589 -0.92158 -4.09914
H 7.64820 0.87464 -3.64958
H 7.71074 -0.54864 -2.58907
H 7.67721 -0.80490 -4.30746
H 12.10162 0.20161 -4.02859
H 12.26326 1.99146 -4.15311
H 10.74652 1.08493 -6.38588
H 11.86985 -0.18805 -6.06571
H 10.42669 -4.30368 -1.42696
H 10.05576 -5.20090 0.00713
H 11.77103 -4.93491 -0.38517
H 10.90584 -5.86455 1.95965
H 12.01917 -5.92036 3.44427
H 12.55574 -5.28496 1.86256
H 10.48205 0.33758 5.70671
H 8.78925 2.73684 5.22222
H 8.51994 -0.28946 6.64803
H 7.83024 -0.34891 4.98767
H 7.68124 1.06640 5.89889
H 10.65124 2.63662 6.92362
H 10.64048 4.07514 5.92085
H 12.36188 1.70509 5.96309
H 12.87048 3.42887 5.84182
H 12.24408 2.51454 4.39754
H 10.21024 7.32058 1.12501
H 9.65238 6.49133 2.61408
H 8.50725 6.66783 1.22491
H 9.22156 4.41385 -3.32635
H 11.53573 5.82268 -6.67936
H 11.98325 4.00028 -6.62895
H 12.82189 5.29039 -5.73198

