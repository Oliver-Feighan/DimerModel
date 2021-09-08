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
Mg 10.50181 0.94860 0.79344
C 10.23195 -1.83828 -1.45145
C 10.54505 2.96727 -1.95569
C 11.19074 3.37983 2.71226
C 11.57318 -1.29818 3.23267
N 10.39515 0.50391 -1.37062
C 10.23197 -0.58743 -2.10321
C 9.97864 -0.31813 -3.61082
C 10.66572 1.12958 -3.71527
C 10.49773 1.58492 -2.26712
C 12.15169 1.06937 -4.02769
C 8.44976 -0.25474 -3.95802
C 8.12728 -0.30348 -5.45314
H 6.91909 -1.00255 -5.88215
N 10.33716 2.90420 0.49834
C 10.35159 3.55261 -0.68515
C 10.21863 4.98174 -0.36255
C 10.29023 5.13756 1.11573
C 10.63725 3.78598 1.50676
C 9.96396 6.03105 -1.40539
C 10.09319 6.31192 2.02581
O 10.43801 6.27292 3.25967
C 9.57112 7.53968 1.40811
N 11.47025 1.06524 2.62942
C 11.59016 2.23890 3.27445
C 12.07243 2.03603 4.69102
C 12.02004 0.45728 4.94074
C 11.58220 -0.00701 3.53429
C 13.47834 2.55281 4.99554
C 11.20527 -0.03441 6.16838
C 9.73362 0.23634 5.95311
N 10.67386 -1.19847 0.96245
C 11.17109 -1.93033 2.06301
C 11.19866 -3.35233 1.66356
C 10.73356 -3.34644 0.33987
C 10.45618 -2.01258 -0.06434
C 11.65934 -4.43330 2.46692
C 10.42424 -4.16417 -0.81092
O 10.38081 -5.37398 -1.03053
C 10.20564 -3.21112 -2.06403
C 8.94498 -3.64297 -2.76172
O 7.84156 -3.53290 -2.19872
O 9.22389 -3.90939 -4.05075
C 8.00264 -4.33129 -4.78721
H 10.79558 3.66996 -2.75089
H 11.34397 4.17015 3.45077
H 11.80148 -2.07136 3.96938
H 10.46781 -1.01135 -4.28456
H 10.11295 1.78164 -4.38656
H 12.52615 0.05663 -4.01112
H 12.63693 1.66469 -3.26534
H 12.35772 1.57040 -4.97788
H 8.06885 0.74583 -3.71617
H 7.91054 -1.02828 -3.44665
H 8.98941 -0.62410 -6.05053
H 7.92562 0.71360 -5.79797
H 10.12596 5.64002 -2.41744
H 10.73704 6.79742 -1.28333
H 8.98277 6.49910 -1.30114
H 10.24331 7.85905 0.61524
H 9.41436 8.23781 2.23268
H 8.61044 7.31120 0.93324
H 11.41363 2.52565 5.39767
H 13.02254 0.03836 5.11203
H 13.50474 3.28177 5.82741
H 13.88392 2.99074 4.09318
H 14.21251 1.78369 5.22926
H 11.49557 0.52190 7.06303
H 11.34444 -1.08827 6.38516
H 9.62829 1.28243 6.25589
H 9.12842 -0.41807 6.57722
H 9.47317 0.17301 4.89018
H 10.94738 -5.23850 2.34339
H 11.75309 -4.14112 3.51076
H 12.63203 -4.61971 2.00388
H 11.08226 -3.30747 -2.70577
H 8.22141 -5.30435 -5.23774
H 7.76680 -3.49871 -5.47638
H 7.11975 -4.56249 -4.20282

