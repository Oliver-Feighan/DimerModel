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
Mg 31.54023 2.83213 2.39832
C 31.74347 0.69630 -0.47847
C 30.76872 0.05399 4.22077
C 30.91292 4.73947 4.85515
C 31.19570 5.51945 0.20622
N 31.31573 0.70597 1.83150
C 31.47738 0.00295 0.72065
C 31.42600 -1.53561 0.91998
C 30.56170 -1.57409 2.27289
C 30.92671 -0.20179 2.83523
C 29.06220 -1.60662 2.02836
C 32.85048 -2.15896 1.13046
C 32.90723 -3.68439 1.02284
H 34.09600 -4.29703 0.43634
N 31.41297 2.42903 4.33737
C 31.11155 1.24167 4.90358
C 31.12599 1.47345 6.35624
C 31.29842 2.93190 6.59685
C 31.19160 3.43691 5.24273
C 31.06392 0.35797 7.35888
C 31.51000 3.73983 7.84132
O 31.39586 5.01647 7.84283
C 31.76271 2.98543 9.07765
N 30.90451 4.80763 2.52328
C 30.75973 5.41322 3.71501
C 30.56338 6.90141 3.55036
C 30.84916 7.20424 2.00628
C 31.08439 5.76245 1.50496
C 29.18147 7.43569 3.92621
C 31.92128 8.28192 1.68714
C 33.28815 7.78981 2.10532
N 31.66107 3.12089 0.26055
C 31.45810 4.32393 -0.45049
C 31.53217 3.99707 -1.88928
C 31.75193 2.61162 -1.91599
C 31.79053 2.10716 -0.58801
C 31.35451 4.91636 -2.96141
C 31.94943 1.45949 -2.76586
O 32.09902 1.28728 -3.97487
C 31.82659 0.14817 -1.87607
C 32.98690 -0.74932 -2.20881
O 34.15095 -0.40283 -1.94136
O 32.51885 -1.95309 -2.58598
C 33.63261 -2.88187 -2.91569
H 30.29871 -0.71111 4.83931
H 30.79835 5.45927 5.66893
H 31.19613 6.31768 -0.53910
H 30.91276 -2.07843 0.13531
H 30.90440 -2.36226 2.93825
H 28.82153 -1.47880 0.98343
H 28.65106 -0.79438 2.61341
H 28.63262 -2.52352 2.44198
H 33.14458 -2.03362 2.18048
H 33.56090 -1.72242 0.45566
H 31.99820 -4.09985 0.57117
H 32.92058 -4.10321 2.03201
H 30.77404 -0.58982 6.88872
H 30.23784 0.58667 8.04091
H 31.98433 0.26101 7.93862
H 30.92748 2.31604 9.26935
H 31.97730 3.73698 9.83988
H 32.64483 2.35302 8.92748
H 31.27269 7.45419 4.15402
H 29.95068 7.57364 1.49042
H 29.21469 8.22660 4.69915
H 28.57301 6.60494 4.25814
H 28.59838 7.83223 3.09678
H 31.72863 9.19010 2.26356
H 31.95078 8.56511 0.64028
H 33.31832 8.02347 3.17377
H 34.06856 8.31995 1.56316
H 33.36242 6.70085 2.00382
H 32.12532 4.69905 -3.68883
H 31.41210 5.94636 -2.61565
H 30.34502 4.64739 -3.28370
H 30.86851 -0.31772 -2.10974
H 33.45672 -3.24300 -3.93369
H 33.64078 -3.63731 -2.10761
H 34.62620 -2.45953 -3.00983

