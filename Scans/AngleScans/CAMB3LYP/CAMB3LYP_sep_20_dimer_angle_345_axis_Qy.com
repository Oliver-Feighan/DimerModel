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
Mg 10.51511 0.95117 0.80682
C 9.67025 -1.34258 3.43437
C 8.72199 3.31171 2.49417
C 10.86113 2.86376 -1.70111
C 11.24344 -1.81469 -1.18448
N 9.38822 0.87387 2.70851
C 9.18873 -0.03472 3.65138
C 8.50049 0.50818 4.93245
C 7.76107 1.77960 4.28777
C 8.69624 2.03361 3.10726
C 6.37930 1.45952 3.74223
C 9.53558 0.94088 6.02943
C 8.93734 1.18904 7.41589
H 9.72872 0.83081 8.58971
N 10.31434 2.91977 0.65562
C 9.56430 3.72371 1.43818
C 9.73892 5.08287 0.90307
C 10.52212 4.99301 -0.35917
C 10.57967 3.55532 -0.53189
C 9.25502 6.30873 1.62177
C 11.10965 6.02978 -1.26807
O 11.54612 5.73206 -2.43584
C 11.07356 7.41889 -0.78789
N 10.78031 0.58205 -1.22191
C 10.95742 1.59255 -2.09108
C 11.40164 1.08001 -3.44029
C 11.72079 -0.47239 -3.22634
C 11.30224 -0.61519 -1.74654
C 10.39308 1.23994 -4.57766
C 13.13289 -0.96492 -3.64592
C 14.18041 -0.36394 -2.73657
N 10.65303 -1.18472 1.10181
C 10.94680 -2.16322 0.12697
C 10.81271 -3.48182 0.77946
C 10.42474 -3.17733 2.09295
C 10.30452 -1.76926 2.24257
C 10.99362 -4.74465 0.14807
C 10.07924 -3.72736 3.38399
O 10.08895 -4.85360 3.87920
C 9.45260 -2.56994 4.27510
C 10.10890 -2.61810 5.62768
O 11.32077 -2.36777 5.75144
O 9.15982 -2.73696 6.57398
C 9.76089 -2.77781 7.93371
H 8.00015 4.05968 2.82337
H 11.09735 3.48270 -2.56977
H 11.54869 -2.72158 -1.71082
H 7.77209 -0.16067 5.37500
H 7.76717 2.62710 4.96826
H 6.16976 0.40058 3.77199
H 6.37612 1.81728 2.72094
H 5.62091 2.04558 4.26922
H 9.90073 1.95000 5.79900
H 10.33367 0.22790 6.10257
H 7.91910 0.79081 7.50276
H 8.81652 2.26624 7.55397
H 8.57136 6.05158 2.44024
H 8.63379 6.86944 0.91499
H 10.07253 6.94941 1.95921
H 10.04291 7.70833 -0.59698
H 11.61888 8.00023 -1.53391
H 11.59810 7.47072 0.17273
H 12.30381 1.58248 -3.76696
H 11.04110 -1.11207 -3.80829
H 10.79135 1.81178 -5.43699
H 9.50563 1.71616 -4.18248
H 9.99664 0.30829 -4.97774
H 13.36769 -0.62442 -4.65740
H 13.23352 -2.04511 -3.63827
H 14.35326 0.62244 -3.17733
H 15.08723 -0.96528 -2.74692
H 13.78283 -0.20873 -1.72688
H 11.56765 -5.35539 0.83218
H 11.49633 -4.64112 -0.81127
H 9.95205 -5.05674 0.03413
H 8.37817 -2.74776 4.33484
H 9.40446 -3.69473 8.41318
H 9.48427 -1.81921 8.41154
H 10.83470 -2.90639 8.00336

