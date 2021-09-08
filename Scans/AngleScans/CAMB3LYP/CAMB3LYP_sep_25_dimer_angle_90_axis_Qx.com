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
Mg 13.13368 1.19028 1.00380
C 15.39114 1.49626 3.77676
C 10.72711 0.30532 3.25321
C 11.24702 0.47671 -1.44544
C 15.92777 0.97581 -1.06811
N 13.17535 0.96606 3.20385
C 14.07286 1.16873 4.15661
C 13.50741 1.09024 5.60011
C 12.22724 0.16606 5.30697
C 11.99982 0.52190 3.83921
C 12.52541 -1.32132 5.39718
C 13.08593 2.49642 6.15423
C 12.81606 2.53979 7.65990
H 13.18019 3.74617 8.39782
N 11.16286 0.97189 0.92184
C 10.33905 0.63143 1.93511
C 8.98750 0.58382 1.35612
C 9.10165 0.76266 -0.11688
C 10.54143 0.72234 -0.27652
C 7.74735 0.46385 2.19340
C 8.08303 0.92810 -1.20372
O 8.39775 0.82918 -2.44227
C 6.68989 1.11644 -0.77340
N 13.52310 0.57313 -0.94295
C 12.52273 0.38254 -1.82094
C 13.05375 0.21142 -3.22420
C 14.61021 0.56936 -3.14119
C 14.73388 0.80958 -1.62054
C 12.88353 -1.17944 -3.83488
C 15.13237 1.66737 -4.10784
C 14.54430 3.00829 -3.73167
N 15.26949 1.40977 1.24510
C 16.26070 1.25302 0.25162
C 17.57156 1.38731 0.91949
C 17.24987 1.59157 2.26977
C 15.83870 1.56491 2.43522
C 18.84233 1.26857 0.28926
C 17.78351 1.81308 3.59449
O 18.90594 2.01439 4.05639
C 16.60778 1.63511 4.64919
C 16.65801 2.79836 5.60146
O 16.42957 3.95336 5.20096
O 16.75183 2.33335 6.86056
C 16.79365 3.44960 7.84242
H 9.96323 -0.20086 3.84409
H 10.63923 0.33453 -2.34201
H 16.84421 1.01897 -1.66052
H 14.15907 0.60610 6.31754
H 11.37487 0.46936 5.90934
H 13.57997 -1.51361 5.52786
H 12.17542 -1.74847 4.46654
H 11.92122 -1.77978 6.18528
H 12.08564 2.74486 5.77683
H 13.81321 3.24136 5.89587
H 13.19441 1.64694 8.17236
H 11.73577 2.50314 7.82001
H 7.98538 0.18400 3.22705
H 7.18043 -0.38906 1.80477
H 7.11961 1.35627 2.14667
H 6.37968 0.26560 -0.17149
H 6.12458 1.30572 -1.68799
H 6.64064 1.99680 -0.12289
H 12.57082 0.89962 -3.90703
H 15.24140 -0.30077 -3.37424
H 12.32580 -1.17133 -4.79036
H 12.38772 -1.81192 -3.11053
H 13.81057 -1.72020 -4.01758
H 14.80402 1.46001 -5.12928
H 16.21417 1.74690 -4.12720
H 13.56472 2.99355 -4.21879
H 15.16264 3.81835 -4.11311
H 14.37399 3.07408 -2.65078
H 19.45854 2.06842 0.67788
H 18.75555 1.32291 -0.79393
H 19.13563 0.27134 0.62804
H 16.76491 0.68294 5.15738
H 17.70005 3.31500 8.44077
H 15.82642 3.41274 8.37803
H 16.94186 4.45134 7.45642

