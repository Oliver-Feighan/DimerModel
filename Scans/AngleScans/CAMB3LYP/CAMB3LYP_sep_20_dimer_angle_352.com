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
Mg 10.51716 0.95195 0.80479
C 10.24343 -0.89778 3.86789
C 9.26399 3.63964 2.49017
C 10.32185 2.47611 -1.97123
C 10.65006 -2.10345 -0.86631
N 9.87511 1.21329 2.90520
C 9.87445 0.45780 3.99313
C 9.53185 1.21687 5.30309
C 8.70549 2.43700 4.66491
C 9.34004 2.46194 3.27588
C 7.22294 2.14268 4.50753
C 10.81370 1.72180 6.05403
C 10.57314 2.19769 7.48841
H 11.60876 1.94847 8.48734
N 10.35615 2.89075 0.41272
C 9.84377 3.84415 1.21881
C 9.93370 5.10220 0.46158
C 10.38903 4.78677 -0.91969
C 10.35255 3.33848 -0.88499
C 9.67932 6.44528 1.08193
C 10.77890 5.64650 -2.08376
O 10.91295 5.16272 -3.26308
C 10.90792 7.08719 -1.82057
N 10.27682 0.29223 -1.15209
C 10.27698 1.16002 -2.17915
C 10.36759 0.43869 -3.50273
C 10.67325 -1.08851 -3.13967
C 10.61548 -0.99878 -1.59889
C 9.12289 0.51204 -4.38688
C 11.92585 -1.73086 -3.79622
C 13.18091 -1.08607 -3.25360
N 10.64559 -1.12817 1.37454
C 10.66306 -2.24847 0.51506
C 10.64195 -3.45327 1.36973
C 10.59001 -2.94631 2.67694
C 10.55910 -1.52607 2.63893
C 10.62184 -4.79987 0.90889
C 10.54358 -3.28984 4.07998
O 10.63127 -4.33599 4.72154
C 10.18944 -1.98093 4.90928
C 11.14768 -1.88995 6.06514
O 12.36213 -1.71020 5.86668
O 10.44849 -1.81229 7.21208
C 11.35510 -1.70926 8.38650
H 8.66871 4.47398 2.86219
H 10.36557 2.95291 -2.95318
H 10.78834 -3.09267 -1.30790
H 8.90696 0.66700 5.99657
H 8.90406 3.36682 5.19157
H 6.98910 1.11489 4.74225
H 6.98863 2.35726 3.47307
H 6.63369 2.84652 5.10227
H 11.14885 2.66216 5.59781
H 11.58032 0.97156 6.04600
H 9.59164 1.88697 7.86663
H 10.52714 3.28943 7.48870
H 9.20226 6.35037 2.06519
H 8.92758 6.94638 0.46253
H 10.57600 7.06740 1.12195
H 9.96365 7.47109 -1.44184
H 11.27956 7.52200 -2.75047
H 11.64822 7.23276 -1.02598
H 11.18288 0.82783 -4.10021
H 9.85195 -1.75257 -3.44679
H 9.32443 0.93228 -5.39031
H 8.37301 1.09838 -3.87267
H 8.60952 -0.43551 -4.54116
H 11.92419 -1.54858 -4.87368
H 11.98695 -2.80429 -3.65129
H 13.27841 -0.18343 -3.86419
H 14.03707 -1.74495 -3.38366
H 13.04178 -0.76737 -2.21402
H 11.32059 -5.35020 1.52484
H 10.88421 -4.86335 -0.14511
H 9.57279 -5.05119 1.08662
H 9.15476 -2.07361 5.24160
H 11.09117 -2.52519 9.06642
H 11.23485 -0.67756 8.76714
H 12.40920 -1.90169 8.22435

