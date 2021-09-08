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
Mg 13.14704 1.18658 0.99291
C 15.71600 -0.98030 2.25154
C 14.91729 3.79524 2.29511
C 10.69142 3.07131 0.29639
C 11.07128 -1.60717 0.81446
N 15.05608 1.26973 2.10687
C 15.97454 0.38845 2.47312
C 17.27682 1.01465 3.04013
C 16.67539 2.42253 3.52494
C 15.49400 2.52800 2.56272
C 16.13375 2.38881 4.94439
C 18.37606 1.20877 1.93740
C 19.77448 1.52546 2.47209
H 20.93029 0.98790 1.75963
N 13.05373 3.16086 0.81384
C 13.86571 4.06977 1.39336
C 13.36787 5.38576 0.96383
C 12.09642 5.18533 0.21688
C 11.88237 3.76892 0.43610
C 14.12557 6.65964 1.20164
C 11.21198 6.11718 -0.55489
O 10.03224 5.77630 -0.92210
C 11.73177 7.47270 -0.78699
N 11.10627 0.83393 0.81058
C 10.26459 1.81728 0.44653
C 8.89731 1.26954 0.11347
C 9.06408 -0.32089 0.09648
C 10.54303 -0.42489 0.52926
C 7.77430 1.65298 1.07690
C 8.61766 -1.06160 -1.19376
C 9.53402 -0.69918 -2.34020
N 13.37991 -0.94408 1.26515
C 12.37497 -1.93133 1.16757
C 12.99092 -3.21876 1.54912
C 14.31605 -2.88482 1.86690
C 14.50675 -1.48481 1.71494
C 12.32224 -4.47358 1.61546
C 15.59407 -3.39673 2.30659
O 16.05695 -4.51824 2.51077
C 16.52347 -2.16782 2.69697
C 17.86802 -2.38066 2.05711
O 17.98766 -2.37038 0.81930
O 18.81927 -2.34394 3.00797
C 20.17169 -2.53925 2.42091
H 15.27426 4.65729 2.85930
H 9.83855 3.65919 -0.05074
H 10.51657 -2.53967 0.69029
H 17.70688 0.48454 3.88144
H 17.37962 3.23266 3.35433
H 16.13532 1.38904 5.35256
H 15.12311 2.77071 4.88275
H 16.68416 3.09313 5.57472
H 18.17105 2.13586 1.38676
H 18.42154 0.35454 1.29023
H 19.85938 1.32681 3.54740
H 19.94423 2.60131 2.38407
H 14.94267 6.51375 1.91890
H 13.44077 7.34948 1.70665
H 14.47353 7.12193 0.27542
H 11.94033 7.94802 0.16853
H 10.99758 7.96096 -1.43077
H 12.68859 7.39479 -1.31523
H 8.57675 1.59976 -0.86699
H 8.47048 -0.80130 0.88811
H 6.92794 2.16339 0.57969
H 8.19103 2.27812 1.85540
H 7.35157 0.82658 1.64573
H 7.61415 -0.74240 -1.48571
H 8.59367 -2.14085 -1.08587
H 9.11992 0.24860 -2.69699
H 9.49822 -1.46220 -3.11516
H 10.55135 -0.50087 -1.98331
H 12.98338 -5.20281 1.16644
H 11.36163 -4.43964 1.10571
H 12.20908 -4.57729 2.69792
H 16.58804 -2.13867 3.78534
H 20.62816 -3.38494 2.94447
H 20.67911 -1.56006 2.50727
H 20.22774 -2.87272 1.39130

