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
Mg 7.89096 0.70999 0.59168
C 7.48490 4.02431 -0.72341
C 7.35239 -0.44962 -2.57058
C 7.79517 -2.22586 1.79129
C 7.26229 2.03691 3.75234
N 7.48322 1.75119 -1.31670
C 7.41299 3.01518 -1.70623
C 7.33744 3.21099 -3.24425
C 6.73156 1.77732 -3.63970
C 7.24661 0.95964 -2.45700
C 5.21248 1.73997 -3.62079
C 8.75057 3.43487 -3.88842
C 8.72622 3.93852 -5.33324
H 9.76384 4.87210 -5.76248
N 8.10098 -1.04977 -0.30146
C 7.86011 -1.33837 -1.59760
C 8.14655 -2.77322 -1.75025
C 8.41752 -3.34343 -0.40260
C 8.08743 -2.20804 0.43526
C 8.22276 -3.44995 -3.08803
C 8.88287 -4.69196 0.05686
O 8.82172 -5.04180 1.28841
C 9.32507 -5.62466 -0.99005
N 7.36739 -0.03269 2.46161
C 7.46497 -1.34531 2.73623
C 7.30135 -1.60697 4.21446
C 7.31464 -0.16684 4.90974
C 7.39772 0.72191 3.64934
C 6.03268 -2.35639 4.62096
C 8.35322 0.05818 6.04260
C 9.75256 0.04072 5.47083
N 7.63410 2.68640 1.42586
C 7.35371 3.01498 2.77029
C 7.15300 4.47743 2.82830
C 7.30905 4.89863 1.49912
C 7.56738 3.77526 0.66796
C 6.82082 5.22599 3.99254
C 7.30345 6.03680 0.60849
O 7.22470 7.25561 0.75693
C 7.29174 5.50608 -0.88971
C 8.33602 6.26668 -1.66009
O 9.54209 6.13968 -1.38466
O 7.75995 6.84640 -2.72894
C 8.75780 7.60353 -3.53074
H 6.97212 -0.92731 -3.47395
H 7.85854 -3.20361 2.27446
H 7.15948 2.53943 4.71648
H 6.66944 3.99950 -3.56953
H 7.15747 1.41161 -4.57041
H 4.79186 2.65395 -3.22853
H 4.94667 0.90060 -2.99162
H 4.82879 1.50309 -4.61728
H 9.23468 2.46106 -4.03704
H 9.34384 4.09469 -3.28568
H 7.73466 4.30519 -5.62528
H 8.90543 3.09098 -5.99925
H 7.81511 -2.81463 -3.88401
H 7.54358 -4.30857 -3.05119
H 9.22834 -3.80465 -3.32396
H 8.51243 -5.78959 -1.69351
H 9.70382 -6.50071 -0.45998
H 10.13957 -5.15674 -1.55419
H 8.12965 -2.19002 4.59799
H 6.35301 0.06147 5.39241
H 6.23672 -3.28626 5.18469
H 5.46186 -2.57098 3.72729
H 5.32528 -1.77914 5.21389
H 8.30438 -0.75713 6.76857
H 8.20391 0.98419 6.58777
H 9.98496 -1.02721 5.41975
H 10.44235 0.56367 6.13019
H 9.76389 0.42992 4.44611
H 7.43792 6.11451 3.97797
H 6.98095 4.64459 4.89813
H 5.75979 5.41376 3.80754
H 6.28944 5.66978 -1.28741
H 8.38627 8.62929 -3.61711
H 8.88211 7.03208 -4.46976
H 9.73390 7.76855 -3.08979

