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
Mg 15.77259 1.42665 1.20213
C 15.37330 0.62403 4.67710
C 14.55474 4.52801 1.93221
C 15.68711 2.01785 -1.91405
C 15.86029 -1.99837 0.56402
N 15.08754 2.33892 3.09706
C 15.03947 1.95837 4.36482
C 14.68670 3.09277 5.36394
C 13.91045 4.07028 4.35359
C 14.57882 3.65094 3.04588
C 12.42448 3.77004 4.24823
C 15.96365 3.78117 5.96186
C 15.70192 4.68306 7.17009
H 16.70572 4.73630 8.22935
N 15.67539 3.15087 0.22427
C 15.17057 4.31691 0.67900
C 15.31391 5.27601 -0.42722
C 15.79348 4.53900 -1.62787
C 15.71567 3.17391 -1.14772
C 15.08230 6.74989 -0.26148
C 16.23526 4.98758 -2.98795
O 16.38417 4.15934 -3.95473
C 16.39813 6.43609 -3.17972
N 15.56123 0.19712 -0.46098
C 15.61048 0.70331 -1.70567
C 15.71282 -0.39464 -2.73739
C 15.96678 -1.73941 -1.90993
C 15.87434 -1.17503 -0.47525
C 14.49238 -0.57543 -3.63974
C 17.21641 -2.57756 -2.29541
C 18.47555 -1.82034 -1.93962
N 15.82898 -0.37615 2.39158
C 15.83586 -1.70793 1.92211
C 15.76042 -2.58750 3.10663
C 15.69111 -1.69913 4.19037
C 15.70087 -0.36043 3.71363
C 15.71379 -4.00999 3.08488
C 15.60120 -1.58949 5.62863
O 15.64408 -2.38648 6.56496
C 15.26389 -0.08138 6.00037
C 16.19609 0.34536 7.10101
O 17.41955 0.43144 6.89518
O 15.47185 0.78835 8.14485
C 16.35235 1.23330 9.25776
H 13.97422 5.44782 2.00866
H 15.76786 2.16558 -2.99347
H 15.98147 -3.07823 0.45494
H 14.03009 2.79721 6.17344
H 14.12217 5.11363 4.57256
H 12.15637 2.87044 4.78206
H 12.22132 3.65755 3.19129
H 11.84095 4.63478 4.57693
H 16.33590 4.52702 5.24779
H 16.70895 3.05094 6.21072
H 14.70326 4.52383 7.59467
H 15.68648 5.72170 6.83096
H 14.57915 6.97382 0.68726
H 14.36005 7.04837 -1.02903
H 15.99480 7.33648 -0.38765
H 15.45610 6.93651 -2.96850
H 16.80414 6.55378 -4.18629
H 17.12281 6.80681 -2.44621
H 16.55287 -0.22571 -3.39987
H 15.13489 -2.45017 -2.02223
H 14.72978 -0.49115 -4.71706
H 13.74699 0.15572 -3.35628
H 13.95658 -1.51410 -3.50928
H 17.24589 -2.73853 -3.37586
H 17.24395 -3.55403 -1.82346
H 18.61298 -1.15370 -2.79623
H 19.31583 -2.50330 -1.83218
H 18.32032 -1.19227 -1.05470
H 16.38177 -4.35533 3.86272
H 15.99966 -4.40231 2.11116
H 14.65412 -4.17368 3.29839
H 14.21930 -0.04660 6.31215
H 16.04935 0.67375 10.14814
H 16.25185 2.33428 9.29632
H 17.40427 0.97995 9.19657

