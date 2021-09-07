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
Mg 31.53905 2.84715 2.39607
C 30.98381 3.83601 5.80091
C 30.44697 5.94457 1.47524
C 31.58803 1.85921 -0.61873
C 31.48432 -0.46016 3.49353
N 30.83440 4.57929 3.57734
C 30.72507 4.86082 4.86697
C 30.39323 6.34554 5.17530
C 29.69950 6.73178 3.77944
C 30.39112 5.71534 2.87326
C 28.20530 6.45601 3.75337
C 31.68060 7.20447 5.43399
C 31.42140 8.58460 6.04211
H 32.38998 9.11797 6.99596
N 31.55784 3.88506 0.70464
C 31.09425 5.13800 0.51358
C 31.32091 5.43834 -0.90860
C 31.80520 4.20058 -1.57802
C 31.64583 3.24055 -0.50423
C 31.15462 6.81384 -1.48632
C 32.31397 3.92388 -2.96034
O 32.45585 2.72785 -3.39876
C 32.55244 5.09421 -3.81763
N 31.32576 0.97241 1.52303
C 31.44151 0.81192 0.19308
C 31.52620 -0.65028 -0.17459
C 31.68718 -1.43299 1.21071
C 31.57311 -0.24291 2.18851
C 30.32991 -1.21383 -0.94120
C 32.90764 -2.38455 1.34371
C 34.18876 -1.58207 1.35653
N 31.46856 1.84386 4.30827
C 31.42763 0.45149 4.53998
C 31.26991 0.25687 5.99589
C 31.20639 1.56004 6.51183
C 31.29645 2.50036 5.45014
C 31.15599 -0.99694 6.66005
C 31.07296 2.35395 7.71208
O 31.04575 2.10852 8.91738
C 30.79579 3.86175 7.29242
C 31.70938 4.74372 8.09864
O 32.94184 4.68830 7.94189
O 30.97209 5.65464 8.75970
C 31.83451 6.55980 9.56514
H 29.90889 6.80104 1.06803
H 31.71244 1.46427 -1.62968
H 31.55729 -1.46073 3.92505
H 29.69607 6.49525 5.99098
H 29.95336 7.74509 3.47925
H 27.87643 5.93403 4.63967
H 28.03306 5.85152 2.87248
H 27.65300 7.38640 3.59315
H 32.11221 7.50209 4.46969
H 32.38118 6.66701 6.04305
H 30.40240 8.67611 6.43733
H 31.46726 9.32961 5.24407
H 30.63079 7.48149 -0.79116
H 30.47403 6.72211 -2.33969
H 32.09793 7.24275 -1.83112
H 31.62881 5.65818 -3.92397
H 32.99763 4.69994 -4.73318
H 33.26866 5.75488 -3.31647
H 32.39543 -0.84435 -0.79094
H 30.82651 -2.08793 1.41106
H 30.60758 -1.66725 -1.91145
H 29.61117 -0.41807 -1.08516
H 29.74565 -1.95828 -0.40302
H 32.96611 -3.04875 0.47793
H 32.87232 -3.01123 2.22858
H 34.38701 -1.41670 0.29331
H 34.99123 -2.14900 1.82415
H 34.03376 -0.60060 1.81949
H 31.77983 -0.93983 7.54213
H 31.45570 -1.81845 6.01270
H 30.08305 -1.00987 6.86961
H 29.74407 4.06960 7.49304
H 31.47495 6.50869 10.59751
H 31.78559 7.54422 9.06287
H 32.87456 6.28176 9.68942

