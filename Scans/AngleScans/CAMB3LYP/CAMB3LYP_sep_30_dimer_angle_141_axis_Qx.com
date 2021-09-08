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
Mg 15.75768 1.42096 1.20367
C 17.07275 3.62472 3.71227
C 15.01691 -0.75923 3.72014
C 15.00663 -0.72436 -1.01015
C 17.60123 3.21910 -1.14446
N 16.05428 1.53196 3.39272
C 16.50790 2.44444 4.23886
C 16.27413 2.10624 5.73571
C 16.16484 0.50924 5.60719
C 15.68120 0.41289 4.16169
C 17.50700 -0.19603 5.70979
C 14.94504 2.73028 6.28899
C 14.80700 2.70224 7.81277
H 14.13781 3.82047 8.47172
N 14.67499 -0.23693 1.33692
C 14.46112 -0.98257 2.44125
C 13.61666 -2.10731 2.00995
C 13.48584 -2.05451 0.52846
C 14.42157 -0.99212 0.21890
C 12.96105 -3.05011 2.97676
C 12.66515 -2.83995 -0.44916
O 12.88629 -2.78588 -1.71051
C 11.65649 -3.74551 0.12010
N 16.39539 1.13254 -0.75401
C 15.87059 0.15637 -1.51533
C 16.27722 0.31396 -2.96108
C 16.98976 1.74270 -3.05369
C 16.94943 2.14272 -1.56246
C 17.21734 -0.75975 -3.50854
C 16.42928 2.74201 -4.10245
C 15.03730 3.17904 -3.70651
N 16.95101 3.22206 1.21134
C 17.65597 3.78284 0.12384
C 18.41161 4.94026 0.64519
C 18.10957 4.95886 2.01520
C 17.24397 3.87654 2.32951
C 19.28019 5.77601 -0.11191
C 18.33464 5.64167 3.26892
O 18.91010 6.67640 3.60345
C 17.77422 4.73337 4.44675
C 16.94905 5.60512 5.35321
O 15.89423 6.12195 4.94475
O 17.42328 5.50884 6.60868
C 16.63042 6.34732 7.54684
H 14.95045 -1.60647 4.40327
H 14.69217 -1.37069 -1.83286
H 18.12204 3.88960 -1.83138
H 17.09248 2.37112 6.39444
H 15.41732 0.10852 6.28689
H 18.32906 0.50381 5.73403
H 17.57449 -0.82878 4.83451
H 17.51346 -0.87073 6.57065
H 14.10308 2.08284 6.01211
H 14.81848 3.73425 5.93310
H 15.75910 2.47782 8.30885
H 14.15837 1.86675 8.08713
H 13.37356 -2.94108 3.98738
H 13.24407 -4.06461 2.67585
H 11.87189 -2.97035 2.96996
H 12.14410 -4.46228 0.77656
H 11.11199 -4.15094 -0.73491
H 10.97365 -3.15973 0.74557
H 15.40955 0.31132 -3.60934
H 18.05154 1.65160 -3.32601
H 16.81584 -1.27833 -4.39958
H 17.42409 -1.46808 -2.71730
H 18.21414 -0.41012 -3.77190
H 16.33662 2.25607 -5.07685
H 17.05192 3.61998 -4.23835
H 14.40699 2.36960 -4.08676
H 14.78601 4.12839 -4.17517
H 14.92625 3.19885 -2.61613
H 19.06938 6.79468 0.18536
H 19.13553 5.63501 -1.18101
H 20.25145 5.40274 0.22364
H 18.63190 4.30129 4.96372
H 17.33469 7.01716 8.04985
H 16.07009 5.63616 8.18248
H 15.93315 7.05725 7.11759

