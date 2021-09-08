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
Mg 7.87457 0.71320 0.60399
C 5.62345 -2.03992 0.12247
C 5.59173 2.66795 -1.00913
C 10.01941 3.05136 0.61123
C 10.40316 -1.62672 1.13006
N 5.89975 0.25772 -0.28195
C 5.11947 -0.81037 -0.35044
C 3.69010 -0.52938 -0.88665
C 4.00420 0.80714 -1.71968
C 5.22044 1.30144 -0.93922
C 4.43354 0.53656 -3.15208
C 2.66144 -0.25882 0.26689
C 1.19279 -0.28295 -0.16237
H 0.20328 -0.81310 0.77160
N 7.64380 2.66842 0.35413
C 6.64703 3.29206 -0.30829
C 6.94490 4.73012 -0.22229
C 8.28525 4.89494 0.40303
C 8.71993 3.51379 0.46269
C 5.96545 5.79011 -0.63525
C 9.05450 6.09779 0.85866
O 10.29963 6.02832 1.15476
C 8.33080 7.37754 0.85847
N 9.95346 0.71947 0.62662
C 10.64077 1.87388 0.68020
C 12.10274 1.62556 0.96545
C 12.21098 0.07304 1.33362
C 10.74264 -0.34529 1.10045
C 13.06708 1.94600 -0.17637
C 12.87431 -0.28351 2.69212
C 11.99999 0.18624 3.83242
N 7.98774 -1.43438 0.81171
C 9.15056 -2.21296 1.00078
C 8.73546 -3.63038 0.96861
C 7.35246 -3.57820 0.73825
C 6.93837 -2.22447 0.61407
C 9.60096 -4.75367 1.09209
C 6.15043 -4.36085 0.56131
O 5.87102 -5.55579 0.64958
C 4.99829 -3.40411 0.02877
C 3.76252 -3.66609 0.84532
O 3.73787 -3.39587 2.05898
O 2.74903 -3.98654 0.02038
C 1.49851 -4.24468 0.78287
H 5.05112 3.31790 -1.69775
H 10.78374 3.82285 0.72993
H 11.11638 -2.41385 1.38387
H 3.29354 -1.29315 -1.54480
H 3.18746 1.52036 -1.64477
H 4.57101 -0.51808 -3.33908
H 5.36612 1.06912 -3.28524
H 3.72470 0.99086 -3.85029
H 2.74710 0.78825 0.58474
H 2.81168 -0.94271 1.07951
H 1.06099 -0.72859 -1.15584
H 0.84958 0.74741 -0.28313
H 5.13271 5.36636 -1.21004
H 6.48216 6.44315 -1.34688
H 5.61572 6.39144 0.20662
H 7.97193 7.58931 -0.14593
H 9.01818 8.10204 1.29946
H 7.44444 7.27831 1.49497
H 12.43641 2.21178 1.81287
H 12.81542 -0.47879 0.59879
H 13.84968 2.67565 0.10507
H 12.49171 2.31073 -1.01688
H 13.57959 1.08556 -0.60320
H 13.82816 0.23937 2.79661
H 13.07338 -1.34338 2.81060
H 12.27274 1.24091 3.93424
H 12.22432 -0.36895 4.74099
H 10.93949 0.14501 3.55804
H 9.10992 -5.45368 1.75494
H 10.57788 -4.46292 1.47269
H 9.64702 -5.07991 0.04964
H 4.84668 -3.63030 -1.02737
H 1.15402 -5.24486 0.50243
H 0.82669 -3.39703 0.55072
H 1.57894 -4.34276 1.85915

