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
Mg 31.52254 2.84180 2.39704
C 29.73675 0.03161 1.05805
C 29.90407 4.77492 0.09961
C 33.40325 5.25180 3.24687
C 33.78674 0.57378 3.76641
N 30.03509 2.35119 0.83529
C 29.38523 1.25975 0.46008
C 28.24159 1.51052 -0.55911
C 28.77730 2.88780 -1.18759
C 29.59712 3.39456 -0.00279
C 29.71517 2.68684 -2.36630
C 26.84954 1.70098 0.13941
C 25.64431 1.64162 -0.80167
H 24.40661 1.04213 -0.31056
N 31.30901 4.79577 2.12325
C 30.59763 5.40869 1.15396
C 30.77576 6.85203 1.37575
C 31.78321 7.03998 2.45483
C 32.22931 5.67366 2.63986
C 29.96845 7.89188 0.65435
C 32.27426 8.25097 3.18871
O 33.32518 8.21406 3.92154
C 31.54256 9.50363 2.95009
N 33.44488 2.92012 3.18496
C 34.00984 4.09507 3.51393
C 35.27449 3.88757 4.31268
C 35.31197 2.32698 4.66017
C 34.05330 1.86652 3.89263
C 36.57497 4.28460 3.61448
C 35.44537 1.94276 6.15925
C 34.19228 2.33807 6.90684
N 31.65127 0.69319 2.58400
C 32.69809 -0.05072 3.17125
C 32.39033 -1.47959 2.95671
C 31.18799 -1.46728 2.23372
C 30.78600 -0.12517 1.99573
C 33.20113 -2.57571 3.36567
C 30.17297 -2.28435 1.60855
O 29.93662 -3.49083 1.56092
C 29.25402 -1.34914 0.71004
C 27.81826 -1.68501 1.00719
O 27.33669 -1.46194 2.13189
O 27.19502 -2.00937 -0.14030
C 25.76525 -2.33988 0.10126
H 29.62479 5.43110 -0.72521
H 34.03364 6.04458 3.65625
H 34.39255 -0.19689 4.24780
H 28.15070 0.75876 -1.33395
H 27.95795 3.56859 -1.40336
H 29.96068 1.64550 -2.51284
H 30.60554 3.25644 -2.13422
H 29.29226 3.14209 -3.26642
H 26.76354 2.73764 0.48967
H 26.72224 0.99253 0.93473
H 25.90779 1.22955 -1.78337
H 25.32195 2.66314 -1.01755
H 29.42597 7.46124 -0.19633
H 30.67955 8.58910 0.19819
H 29.30618 8.44829 1.32103
H 31.56854 9.74055 1.88912
H 31.98521 10.23460 3.62953
H 30.48995 9.34937 3.21243
H 35.24566 4.45264 5.23625
H 36.16917 1.82489 4.18800
H 37.16448 5.03009 4.18088
H 36.33251 4.66047 2.62935
H 37.24795 3.45944 3.38773
H 36.26860 2.49444 6.61979
H 35.63606 0.88680 6.31918
H 34.35910 3.39705 7.12552
H 34.09254 1.75711 7.82149
H 33.31002 2.27009 6.25981
H 32.53399 -3.31709 3.78492
H 33.95511 -2.26548 4.08612
H 33.64227 -2.86046 2.40672
H 29.51192 -1.54039 -0.33231
H 25.59491 -3.34020 -0.30866
H 25.18711 -1.50814 -0.34338
H 25.44893 -2.47573 1.12885

