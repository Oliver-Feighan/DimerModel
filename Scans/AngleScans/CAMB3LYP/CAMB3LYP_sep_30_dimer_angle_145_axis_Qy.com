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
Mg 15.76474 1.42177 1.18745
C 17.75366 -1.10803 -0.40107
C 17.74970 3.73371 -0.34534
C 14.32574 3.65845 2.91775
C 14.70638 -1.01974 3.43792
N 17.49557 1.21597 -0.17422
C 18.12324 0.20576 -0.75718
C 19.15341 0.62696 -1.83926
C 19.45956 2.11337 -1.31415
C 18.14241 2.39426 -0.59382
C 20.58395 2.17141 -0.29350
C 18.53347 0.64325 -3.28058
C 19.55246 0.74366 -4.41783
H 19.27896 0.03047 -5.66243
N 15.66513 3.38582 0.92075
C 16.54436 4.16198 0.25305
C 16.03185 5.53655 0.36285
C 14.86699 5.53027 1.28901
C 14.92149 4.16500 1.77187
C 16.59193 6.67977 -0.43262
C 13.87660 6.58296 1.68546
O 13.09128 6.42588 2.68618
C 13.91739 7.84088 0.92557
N 14.84006 1.38677 3.04919
C 14.24738 2.48602 3.54749
C 13.41557 2.15365 4.76323
C 13.37876 0.55644 4.83503
C 14.31045 0.22821 3.64777
C 13.92672 2.71089 6.09161
C 11.97834 -0.11380 4.88290
C 11.26025 0.09671 3.56934
N 15.99351 -0.71778 1.38252
C 15.48936 -1.53893 2.41481
C 16.00130 -2.90374 2.17474
C 16.79275 -2.77621 1.02329
C 16.79160 -1.42316 0.58882
C 15.75755 -4.04016 2.99647
C 17.63749 -3.47647 0.08268
O 17.93747 -4.65476 -0.10543
C 18.39934 -2.40329 -0.80854
C 18.27818 -2.82431 -2.24749
O 17.16934 -2.83948 -2.81038
O 19.50892 -2.92931 -2.78099
C 19.44066 -3.33393 -4.21048
H 18.44796 4.53360 -0.59334
H 13.72457 4.36878 3.49004
H 14.34283 -1.85630 4.03847
H 20.06500 0.04161 -1.85089
H 19.59636 2.80435 -2.14191
H 20.91353 1.18563 -0.00066
H 20.18113 2.70413 0.55796
H 21.40559 2.78601 -0.67234
H 17.99406 1.58778 -3.42769
H 17.90366 -0.21138 -3.43412
H 20.57422 0.54029 -4.07491
H 19.58673 1.78000 -4.76256
H 17.54739 6.41371 -0.90135
H 16.85011 7.46969 0.28098
H 15.87875 7.07777 -1.15763
H 14.90669 8.28371 1.01303
H 13.07778 8.43284 1.29503
H 13.76705 7.61354 -0.13572
H 12.40375 2.52386 4.65251
H 13.87699 0.17764 5.73953
H 13.18420 3.34285 6.61446
H 14.83401 3.26868 5.90126
H 14.26154 1.96268 6.80804
H 11.36044 0.35458 5.65295
H 12.01503 -1.17585 5.10131
H 10.82690 1.09498 3.68222
H 10.48677 -0.65598 3.43069
H 11.97004 0.13789 2.73498
H 15.54093 -4.86503 2.33077
H 14.93776 -3.86145 3.68914
H 16.71930 -4.11918 3.51034
H 19.43756 -2.37838 -0.47530
H 20.04980 -4.23694 -4.31636
H 19.75430 -2.44415 -4.78812
H 18.48516 -3.67779 -4.58896

