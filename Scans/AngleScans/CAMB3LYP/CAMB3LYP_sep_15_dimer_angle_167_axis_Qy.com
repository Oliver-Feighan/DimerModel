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
Mg 7.87802 0.71191 0.58981
C 9.05078 -1.94245 -1.52160
C 9.23748 2.89427 -1.65157
C 7.25806 3.05110 2.64195
C 7.63938 -1.62699 3.16242
N 8.97624 0.39551 -1.30388
C 9.30880 -0.65554 -2.03806
C 9.88134 -0.30693 -3.43814
C 10.41151 1.17493 -3.11897
C 9.46362 1.53279 -1.97612
C 11.83436 1.20327 -2.58602
C 8.77449 -0.29388 -4.55022
C 9.30504 -0.26635 -5.98532
H 8.56695 -0.99356 -7.01433
N 7.75602 2.67205 0.30429
C 8.35364 3.39138 -0.66870
C 7.96631 4.79042 -0.43031
C 7.22572 4.85898 0.85870
C 7.40655 3.50374 1.33902
C 8.23324 5.88835 -1.41860
C 6.48895 5.96520 1.55129
O 6.12310 5.86716 2.77578
C 6.29069 7.20264 0.78257
N 7.70458 0.76143 2.66100
C 7.37642 1.89789 3.30029
C 7.04055 1.63174 4.74830
C 6.97681 0.03997 4.88899
C 7.39280 -0.35761 3.45571
C 8.02481 2.19402 5.77371
C 5.67021 -0.56312 5.47356
C 4.52613 -0.34871 4.50893
N 8.08742 -1.43110 0.76795
C 7.97109 -2.20471 1.94357
C 8.31018 -3.59688 1.58396
C 8.62509 -3.53198 0.21827
C 8.51124 -2.19035 -0.23623
C 8.34699 -4.70215 2.48014
C 9.03824 -4.29156 -0.93976
O 9.20617 -5.48643 -1.18023
C 9.45480 -3.27511 -2.08854
C 8.79678 -3.72195 -3.36519
O 7.55871 -3.69802 -3.48012
O 9.73948 -3.89597 -4.30934
C 9.13465 -4.32869 -5.59726
H 9.82263 3.65506 -2.16893
H 6.93566 3.80120 3.36781
H 7.49387 -2.43216 3.88575
H 10.70324 -0.93420 -3.76186
H 10.25751 1.84021 -3.96458
H 12.21392 0.20998 -2.39774
H 11.79293 1.77298 -1.66691
H 12.47936 1.77037 -3.26330
H 8.25233 0.67122 -4.52440
H 8.10289 -1.12148 -4.42907
H 10.37343 -0.50938 -6.03453
H 9.24607 0.75934 -6.35747
H 8.93849 5.56774 -2.19528
H 8.76395 6.68108 -0.88028
H 7.31737 6.30288 -1.84523
H 7.25738 7.60076 0.48346
H 5.66786 7.84111 1.41210
H 5.75165 6.95897 -0.13975
H 6.07287 2.04609 5.00331
H 7.75991 -0.34138 5.56068
H 7.55014 2.87133 6.50856
H 8.81689 2.70460 5.24224
H 8.57379 1.44716 6.34469
H 5.39672 -0.04947 6.39850
H 5.74759 -1.62060 5.70322
H 4.20026 0.67092 4.73533
H 3.73015 -1.06747 4.69287
H 4.87911 -0.35919 3.47117
H 7.87133 -5.53066 1.97237
H 7.84729 -4.47010 3.41834
H 9.42706 -4.81439 2.60713
H 10.54292 -3.29114 -2.16143
H 9.62972 -5.26121 -5.88538
H 9.24413 -3.46756 -6.28298
H 8.09550 -4.63604 -5.58475

