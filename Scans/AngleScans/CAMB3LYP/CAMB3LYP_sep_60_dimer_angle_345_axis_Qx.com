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
Mg 31.53972 2.84185 2.40252
C 30.69731 0.69411 5.15145
C 30.87614 5.48275 4.45676
C 31.83693 4.74815 -0.11644
C 31.02857 0.11369 0.29501
N 30.88027 2.94702 4.51115
C 30.65169 2.06649 5.47384
C 30.43658 2.69076 6.87868
C 29.95048 4.15126 6.42095
C 30.63688 4.22143 5.05838
C 28.44922 4.24562 6.20492
C 31.76490 2.76208 7.71083
C 31.57657 3.07222 9.19747
H 32.47716 2.44277 10.15934
N 31.85701 4.79837 2.30294
C 31.54610 5.72172 3.23677
C 31.96406 7.01644 2.67717
C 32.39615 6.80148 1.26940
C 32.01891 5.41424 1.08679
C 32.00379 8.27884 3.48847
C 33.03065 7.69989 0.25133
O 33.10188 7.37128 -0.98551
C 33.48212 9.01523 0.72835
N 31.24221 2.54464 0.36647
C 31.49387 3.52448 -0.51925
C 31.47381 2.99867 -1.93471
C 31.39499 1.40700 -1.80068
C 31.28848 1.28737 -0.26453
C 30.32358 3.49286 -2.81171
C 32.49073 0.58404 -2.53195
C 33.83557 0.82579 -1.88516
N 31.14935 0.73205 2.65629
C 30.94303 -0.22561 1.63924
C 30.59987 -1.49954 2.30390
C 30.60787 -1.18654 3.67155
C 30.91324 0.19035 3.84595
C 30.28541 -2.72239 1.64674
C 30.41698 -1.70537 5.00687
O 30.22661 -2.82258 5.48576
C 30.34213 -0.48173 6.01855
C 31.24023 -0.78831 7.18557
O 32.47080 -0.87785 7.02968
O 30.52641 -0.70997 8.32341
C 31.37635 -0.99527 9.51001
H 30.47691 6.37345 4.94272
H 32.03695 5.33561 -1.01552
H 30.95095 -0.80653 -0.28794
H 29.66851 2.21324 7.47519
H 30.33725 4.92174 7.08286
H 27.96986 3.28094 6.28075
H 28.32013 4.65521 5.21163
H 28.01484 4.97779 6.89163
H 32.33236 3.65012 7.40405
H 32.33335 1.85929 7.59922
H 30.53308 2.95465 9.51409
H 31.78834 4.13175 9.36046
H 31.47269 8.16162 4.44114
H 31.42037 9.02821 2.94265
H 33.01867 8.65618 3.63098
H 32.63793 9.55582 1.14973
H 33.98824 9.47650 -0.12192
H 34.19822 8.86471 1.54393
H 32.38452 3.26257 -2.45823
H 30.45585 1.00942 -2.21271
H 30.66453 3.99131 -3.73873
H 29.70938 4.16208 -2.22395
H 29.60972 2.72757 -3.11158
H 32.58194 0.91309 -3.57001
H 32.29118 -0.48222 -2.54845
H 34.17060 1.75645 -2.35283
H 34.51808 0.00673 -2.10276
H 33.72711 1.01698 -0.81122
H 30.80597 -3.50589 2.18120
H 30.57414 -2.69587 0.59807
H 29.19959 -2.73819 1.77326
H 29.30036 -0.37063 6.32179
H 30.89688 -1.81272 10.05734
H 31.48428 -0.02986 10.03925
H 32.36211 -1.40910 9.33260

