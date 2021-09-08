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
Mg 13.12840 1.18547 0.99664
C 11.95203 -1.63510 -0.88474
C 12.22541 3.13826 -1.65022
C 14.56129 3.62410 2.43446
C 14.94438 -1.05390 2.95448
N 12.21558 0.69818 -0.95824
C 11.75553 -0.40041 -1.53772
C 10.97362 -0.15324 -2.85571
C 11.62406 1.25824 -3.26009
C 12.01811 1.75303 -1.86991
C 12.88694 1.12199 -4.09430
C 9.42772 -0.02649 -2.61879
C 8.57503 -0.09343 -3.88773
H 7.27007 -0.74487 -3.81585
N 12.93522 3.13956 0.70894
C 12.53486 3.76010 -0.42063
C 12.58047 5.20070 -0.12607
C 13.19767 5.38588 1.21538
C 13.61682 4.02865 1.50233
C 11.99610 6.23726 -1.04121
C 13.39150 6.58810 2.08904
O 14.16549 6.56056 3.11032
C 12.72182 7.82477 1.66048
N 14.70918 1.29774 2.34233
C 15.09962 2.47861 2.85322
C 16.06292 2.28509 3.99995
C 16.05099 0.71643 4.31101
C 15.10929 0.24184 3.18258
C 17.49922 2.74252 3.74663
C 15.72975 0.29065 5.76984
C 14.29293 0.62475 6.10047
N 13.27525 -0.96298 1.17170
C 14.11734 -1.69252 2.03930
C 13.94573 -3.12280 1.71179
C 13.02573 -3.12472 0.65251
C 12.66579 -1.78873 0.32839
C 14.63202 -4.20590 2.32988
C 12.28530 -3.95258 -0.27221
O 12.12153 -5.16364 -0.41453
C 11.65350 -3.01851 -1.39235
C 10.21008 -3.40664 -1.56096
O 9.39670 -3.23282 -0.63641
O 9.98437 -3.71437 -2.85121
C 8.56360 -4.09525 -3.07071
H 12.18948 3.81066 -2.50781
H 15.00372 4.42267 3.03448
H 15.40104 -1.82028 3.58440
H 11.15518 -0.88329 -3.63534
H 10.88591 1.92022 -3.70537
H 13.20535 0.09354 -4.17822
H 13.63971 1.71089 -3.58677
H 12.74543 1.59162 -5.07193
H 9.19820 0.99586 -2.29191
H 9.08834 -0.76252 -1.91614
H 9.14423 -0.46705 -4.74759
H 8.29621 0.92401 -4.17252
H 11.75959 5.81668 -2.02639
H 12.78606 6.96934 -1.24094
H 11.13960 6.75267 -0.60141
H 13.06494 8.09478 0.66458
H 12.90482 8.54762 2.45790
H 11.64635 7.63075 1.58098
H 15.72881 2.82045 4.88022
H 17.03063 0.25523 4.11755
H 17.85613 3.48779 4.48234
H 17.55849 3.14090 2.74248
H 18.24038 1.94544 3.72306
H 16.34887 0.85262 6.47342
H 15.90201 -0.76346 5.96001
H 14.34345 1.68109 6.38096
H 13.93803 0.01327 6.92762
H 13.65675 0.54997 5.21090
H 13.89693 -4.97964 2.50725
H 15.11435 -3.89523 3.25433
H 15.35818 -4.44763 1.54915
H 12.22767 -3.16982 -2.30724
H 8.56651 -5.08723 -3.53294
H 8.11957 -3.26814 -3.65592
H 7.95098 -4.27205 -2.19442

