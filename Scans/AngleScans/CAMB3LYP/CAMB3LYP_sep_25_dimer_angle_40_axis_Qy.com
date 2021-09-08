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
Mg 13.14811 1.18734 0.99899
C 14.65785 -0.93150 3.47088
C 13.60275 3.78419 3.16331
C 11.42335 2.99181 -0.95971
C 11.80388 -1.68673 -0.44269
N 14.02874 1.28771 3.02549
C 14.62614 0.42978 3.83876
C 15.29390 1.07900 5.08061
C 14.41128 2.41849 5.15552
C 14.01344 2.53160 3.68528
C 13.14585 2.25997 5.98191
C 16.81062 1.40305 4.84279
C 17.59981 1.75087 6.10689
H 18.99264 1.32172 6.19800
N 13.01840 3.16223 0.85222
C 13.26056 4.06800 1.82290
C 13.00437 5.38426 1.21794
C 12.43251 5.17324 -0.13969
C 12.24579 3.73621 -0.12663
C 13.37342 6.67460 1.89043
C 12.09595 6.11127 -1.25908
O 11.38700 5.73859 -2.25975
C 12.54173 7.50365 -1.10406
N 11.63738 0.74754 -0.35991
C 11.09237 1.71081 -1.12345
C 10.22886 1.12039 -2.21255
C 10.50167 -0.45486 -2.17057
C 11.44944 -0.51517 -0.95259
C 8.72794 1.38025 -2.08577
C 10.95700 -1.12580 -3.49540
C 12.33824 -0.64010 -3.87176
N 13.34844 -0.94289 1.29904
C 12.67473 -1.96958 0.60182
C 13.05196 -3.24822 1.23837
C 13.90769 -2.87185 2.28459
C 14.03750 -1.45696 2.31166
C 12.57506 -4.53619 0.86447
C 14.72227 -3.34902 3.37887
O 15.06654 -4.45767 3.78645
C 15.14408 -2.10553 4.27453
C 16.62042 -2.20681 4.54416
O 17.44100 -2.10509 3.61523
O 16.82782 -2.18884 5.87345
C 18.27817 -2.27592 6.19081
H 13.49068 4.62180 3.85232
H 10.89134 3.55933 -1.72676
H 11.50446 -2.63473 -0.89470
H 15.19083 0.51395 5.99924
H 15.01381 3.27235 5.45426
H 12.98876 1.23569 6.28574
H 12.33601 2.59472 5.34697
H 13.16380 2.94542 6.83406
H 16.89293 2.35470 4.30210
H 17.29525 0.59889 4.32405
H 17.05457 1.48283 7.02002
H 17.70145 2.83732 6.16458
H 13.62416 6.52012 2.94715
H 12.46967 7.29314 1.91333
H 14.15918 7.21647 1.35975
H 12.11198 7.92154 -0.19678
H 12.28719 7.99829 -2.04343
H 13.62909 7.50991 -0.96910
H 10.51786 1.50131 -3.18440
H 9.59806 -1.01763 -1.89378
H 8.29531 1.88128 -2.97227
H 8.55768 1.96995 -1.19482
H 8.12055 0.49660 -1.89751
H 10.29285 -0.83726 -4.31381
H 10.96156 -2.20979 -3.45127
H 12.13672 0.30842 -4.37851
H 12.82480 -1.34856 -4.53913
H 12.93375 -0.41680 -2.97892
H 13.43012 -5.19901 0.87180
H 12.09608 -4.51470 -0.11220
H 11.85823 -4.71994 1.66925
H 14.55641 -2.14852 5.19234
H 18.40797 -3.13214 6.85993
H 18.55783 -1.28112 6.58568
H 18.95328 -2.53425 5.38338

