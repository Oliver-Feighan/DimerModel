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
Mg 13.14667 1.17759 0.99368
C 13.07539 2.48354 -2.34824
C 12.64463 -1.91160 -0.36260
C 12.75490 0.05839 3.93670
C 12.50072 4.41092 2.12239
N 12.90752 0.50130 -1.09861
C 12.95994 1.07903 -2.28934
C 12.96513 0.08902 -3.48492
C 12.25798 -1.16049 -2.76553
C 12.65292 -0.86583 -1.31989
C 10.74167 -1.13127 -2.86186
C 14.41764 -0.26782 -3.95916
C 14.49578 -0.98133 -5.31066
H 15.62508 -0.68746 -6.18862
N 13.24516 -0.68012 1.68491
C 13.03503 -1.81861 0.99139
C 13.20607 -2.91525 1.95693
C 13.37192 -2.32480 3.31287
C 13.10377 -0.92956 3.02741
C 13.28076 -4.35708 1.54587
C 13.70288 -2.91956 4.64813
O 13.56140 -2.25504 5.73497
C 14.10837 -4.33273 4.65400
N 12.48527 2.05410 2.75908
C 12.46078 1.35842 3.90945
C 12.21444 2.26624 5.09078
C 12.31967 3.75548 4.51747
C 12.52943 3.43630 3.02097
C 10.87152 2.08737 5.79857
C 13.32545 4.70976 5.21776
C 14.74095 4.24189 4.96708
N 13.02220 3.14097 0.10054
C 12.71472 4.35603 0.75105
C 12.63527 5.39992 -0.29127
C 12.88128 4.71000 -1.48792
C 13.07873 3.32979 -1.21308
C 12.31899 6.77056 -0.07349
C 13.00832 4.83366 -2.92224
O 13.02608 5.77450 -3.71485
C 13.01441 3.37625 -3.55653
C 14.15000 3.30252 -4.54016
O 15.32865 3.38515 -4.15213
O 13.66960 2.93074 -5.74073
C 14.75998 2.83290 -6.74737
H 12.26349 -2.88774 -0.66387
H 12.71580 -0.25570 4.98231
H 12.40047 5.46142 2.40359
H 12.37981 0.40601 -4.33972
H 12.69002 -2.10203 -3.09448
H 10.38332 -0.21045 -3.29755
H 10.38038 -1.23543 -1.84713
H 10.38122 -2.01145 -3.40194
H 14.82403 -1.05281 -3.30850
H 15.03863 0.60667 -3.97578
H 13.55151 -0.91865 -5.86507
H 14.63073 -2.05072 -5.13120
H 12.96136 -4.49585 0.50561
H 12.53095 -4.89666 2.13451
H 14.26195 -4.79832 1.73362
H 13.31471 -4.93796 4.22242
H 14.38985 -4.55170 5.68584
H 14.98204 -4.44872 4.00295
H 12.97418 2.12752 5.85016
H 11.36136 4.29000 4.59347
H 10.97313 1.86106 6.87677
H 10.32235 1.30264 5.29542
H 10.19103 2.93338 5.71858
H 13.17818 4.68698 6.30033
H 13.23188 5.74337 4.90170
H 14.88478 3.47092 5.72998
H 15.44419 5.06199 5.09721
H 14.82753 3.75686 3.98785
H 13.00858 7.34848 -0.67430
H 12.39193 7.03292 0.97989
H 11.28608 6.79159 -0.43118
H 12.04680 3.22367 -4.03616
H 14.47981 3.47831 -7.58560
H 14.87482 1.75327 -6.95940
H 15.72730 3.24172 -6.47979

