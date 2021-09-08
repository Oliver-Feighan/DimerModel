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
Mg 13.14820 1.18714 0.99711
C 15.05371 -0.93159 3.17864
C 14.05101 3.80187 2.99279
C 11.12176 3.00323 -0.63469
C 11.50202 -1.67531 -0.11741
N 14.39522 1.29492 2.82073
C 15.11792 0.43485 3.52239
C 16.01926 1.08556 4.60577
C 15.19244 2.44522 4.82155
C 14.52805 2.54753 3.45013
C 14.10174 2.32590 5.87296
C 17.47039 1.37262 4.08248
C 18.48934 1.71964 5.17018
H 19.86597 1.26089 5.00616
N 13.03146 3.16230 0.84315
C 13.46893 4.07531 1.73551
C 13.12927 5.38876 1.16661
C 12.30885 5.17242 -0.05602
C 12.10015 3.74020 0.01660
C 13.64283 6.67944 1.73568
C 11.78641 6.10269 -1.10857
O 10.89522 5.73256 -1.95214
C 12.28020 7.48677 -1.06386
N 11.40110 0.76295 -0.04694
C 10.74123 1.72783 -0.71139
C 9.67745 1.14228 -1.60910
C 9.92277 -0.43794 -1.59177
C 11.08091 -0.50307 -0.57203
C 8.23228 1.43710 -1.20791
C 10.10837 -1.13645 -2.96666
C 11.40358 -0.68663 -3.60342
N 13.36005 -0.94283 1.29106
C 12.54776 -1.96349 0.75020
C 13.01295 -3.24160 1.32681
C 14.05690 -2.87043 2.18756
C 14.21682 -1.45854 2.16539
C 12.44954 -4.52352 1.07116
C 15.05297 -3.35095 3.11795
O 15.44610 -4.46145 3.47290
C 15.65931 -2.10532 3.89713
C 17.15778 -2.23578 3.88714
O 17.79129 -2.16472 2.81933
O 17.61124 -2.20471 5.15354
C 19.09342 -2.31975 5.19505
H 14.08643 4.65088 3.67602
H 10.46634 3.57215 -1.29819
H 11.10483 -2.62237 -0.48889
H 16.07949 0.53524 5.53706
H 15.85672 3.28938 4.98735
H 13.98467 1.30951 6.21848
H 13.19376 2.67007 5.39529
H 14.29251 3.02209 6.69471
H 17.46817 2.31491 3.51960
H 17.83343 0.55104 3.49611
H 18.12005 1.47598 6.17377
H 18.62098 2.80423 5.18904
H 14.08437 6.53352 2.72924
H 12.77156 7.31816 1.91678
H 14.32538 7.19658 1.05793
H 12.03649 7.92616 -0.09950
H 11.86352 7.97437 -1.94733
H 13.37349 7.47065 -1.13511
H 9.78626 1.50368 -2.62430
H 9.07644 -0.97675 -1.14094
H 7.65076 1.93576 -2.00619
H 8.24366 2.04227 -0.31115
H 7.65403 0.56976 -0.89398
H 9.30814 -0.84416 -3.65103
H 10.10017 -2.21959 -2.90545
H 11.12892 0.25929 -4.07975
H 11.74247 -1.41457 -4.33779
H 12.16025 -0.46474 -2.84193
H 13.27781 -5.20503 0.92963
H 11.79631 -4.50441 0.20127
H 11.89304 -4.68055 1.99902
H 15.25358 -2.12299 4.90938
H 19.32989 -3.16965 5.84270
H 19.46139 -1.32624 5.51331
H 19.59990 -2.60377 4.28003

