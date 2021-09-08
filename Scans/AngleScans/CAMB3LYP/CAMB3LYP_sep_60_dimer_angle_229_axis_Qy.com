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
Mg 31.52233 2.84157 2.39517
C 30.08189 0.02093 0.70752
C 30.32031 4.78243 -0.13917
C 33.15651 5.27246 3.61487
C 33.53978 0.59445 4.13472
N 30.36504 2.34934 0.57565
C 29.82139 1.25217 0.07065
C 28.88363 1.49835 -1.14161
C 29.49684 2.89638 -1.63998
C 30.06874 3.39836 -0.31574
C 30.64314 2.73397 -2.62434
C 27.38149 1.64978 -0.71422
C 26.37550 1.58338 -1.86546
H 25.08136 0.95169 -1.62335
N 31.32049 4.79549 2.11226
C 30.78982 5.41118 1.03495
C 30.89113 6.85366 1.30545
C 31.67418 7.04239 2.55682
C 32.10799 5.68195 2.80403
C 30.21034 7.89000 0.45932
C 31.99202 8.24979 3.38589
O 32.88762 8.22056 4.30235
C 31.29031 9.49201 3.03094
N 33.26082 2.94343 3.53078
C 33.72791 4.12323 3.97556
C 34.82483 3.92645 4.99459
C 34.83124 2.36112 5.32209
C 33.74917 1.88973 4.32593
C 36.22393 4.36055 4.55817
C 34.68997 1.95373 6.81432
C 33.31055 2.31176 7.31870
N 31.66149 0.69295 2.57430
C 32.59605 -0.04061 3.33757
C 32.36576 -1.47125 3.05006
C 31.32016 -1.46971 2.11453
C 30.94014 -0.13171 1.82325
C 33.10976 -2.55838 3.58924
C 30.45864 -2.29533 1.29916
O 30.26227 -3.50515 1.19195
C 29.70371 -1.36269 0.25673
C 28.24553 -1.73134 0.27469
O 27.55697 -1.53707 1.29191
O 27.85564 -2.04782 -0.97359
C 26.41364 -2.41003 -1.00904
H 30.18597 5.44722 -0.99291
H 33.68128 6.07009 4.14580
H 34.06167 -0.17256 4.71091
H 28.95623 0.75847 -1.92969
H 28.71749 3.56482 -1.99656
H 30.93484 1.70025 -2.73607
H 31.46139 3.31658 -2.22174
H 30.38631 3.19643 -3.58167
H 27.20838 2.67838 -0.37253
H 27.12324 0.92537 0.03356
H 26.82731 1.19349 -2.78568
H 26.07667 2.60205 -2.12434
H 29.84652 7.46371 -0.48371
H 30.97863 8.60861 0.15403
H 29.42270 8.42191 0.99725
H 31.50930 9.74767 1.99694
H 31.58150 10.21955 3.79105
H 30.21089 9.31291 3.08902
H 34.61095 4.47483 5.90379
H 35.77266 1.88392 5.01251
H 36.68018 5.10740 5.23502
H 36.16200 4.74862 3.55015
H 36.94564 3.55259 4.45075
H 35.39989 2.51319 7.42847
H 34.87075 0.89904 6.99310
H 33.40987 3.36984 7.57890
H 33.05418 1.71328 8.19056
H 32.56684 2.23792 6.51677
H 32.39254 -3.31965 3.86594
H 33.70835 -2.24613 4.44245
H 33.72894 -2.81795 2.72639
H 30.15649 -1.53089 -0.72116
H 26.34538 -3.40623 -1.45694
H 25.91068 -1.58207 -1.54313
H 25.91353 -2.56969 -0.06096

