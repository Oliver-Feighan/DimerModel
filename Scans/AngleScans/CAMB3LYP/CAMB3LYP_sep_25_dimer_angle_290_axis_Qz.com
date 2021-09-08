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
Mg 13.14522 1.19100 0.99685
C 12.54349 3.00734 4.03295
C 12.13959 3.97765 -0.69367
C 13.23311 -0.51893 -1.67445
C 12.99643 -1.73631 2.88218
N 12.45495 3.17431 1.69130
C 12.32593 3.77059 2.86695
C 12.02041 5.29041 2.78721
C 11.36339 5.32708 1.32234
C 12.05073 4.10560 0.71552
C 9.86437 5.07762 1.32868
C 13.32075 6.16588 2.85560
C 13.07955 7.65811 3.09396
H 14.04020 8.39705 3.90856
N 13.22055 1.77336 -0.89886
C 12.78830 2.94623 -1.40760
C 13.04992 2.87834 -2.85367
C 13.52045 1.50502 -3.18126
C 13.31869 0.84614 -1.90640
C 12.92527 4.06856 -3.75990
C 14.05068 0.88385 -4.43792
O 14.17523 -0.38589 -4.56078
C 14.33174 1.79911 -5.55362
N 12.90868 -0.83865 0.61354
C 13.04755 -1.32791 -0.63116
C 13.10771 -2.83677 -0.62062
C 13.22378 -3.25135 0.91959
C 13.11609 -1.85318 1.56682
C 11.91503 -3.55451 -1.25209
C 14.42030 -4.15910 1.31582
C 15.71809 -3.39957 1.16009
N 13.01455 0.69811 3.09639
C 12.93865 -0.59150 3.66669
C 12.74756 -0.41393 5.12078
C 12.70213 0.97763 5.29393
C 12.83393 1.62156 4.03398
C 12.59304 -1.46024 6.07340
C 12.56198 2.04807 6.25477
O 12.50526 2.11178 7.48215
C 12.32621 3.40763 5.46584
C 13.24245 4.44812 6.04928
O 14.47631 4.33548 5.94197
O 12.51208 5.50697 6.44394
C 13.37770 6.57053 7.01958
H 11.62858 4.71389 -1.31472
H 13.36913 -1.15569 -2.55167
H 13.03891 -2.59852 3.55111
H 11.31046 5.65021 3.52225
H 11.64525 6.22914 0.78549
H 9.50649 4.79880 2.30864
H 9.69669 4.27523 0.62222
H 9.33580 5.94733 0.92798
H 13.77799 6.20637 1.85859
H 13.99710 5.78628 3.59661
H 12.05506 7.86178 3.42845
H 13.15762 8.17948 2.13683
H 12.40228 4.89695 -3.26629
H 12.26008 3.77767 -4.58013
H 13.88442 4.38256 -4.17714
H 13.42294 2.33347 -5.82007
H 14.78650 1.18165 -6.33067
H 15.05200 2.55235 -5.21527
H 13.98468 -3.19252 -1.14735
H 12.34520 -3.82159 1.25538
H 12.20216 -4.24020 -2.07152
H 11.21684 -2.80842 -1.60764
H 11.30402 -4.13153 -0.56006
H 14.48163 -5.01926 0.64462
H 14.35359 -4.54440 2.32780
H 15.94118 -3.50807 0.09447
H 16.49847 -3.84466 1.77404
H 15.57526 -2.33123 1.35986
H 13.20017 -1.19479 6.92865
H 12.88771 -2.42206 5.65886
H 11.51610 -1.40317 6.25283
H 11.27545 3.67588 5.58213
H 12.99648 6.78455 8.02285
H 13.36032 7.39905 6.28677
H 14.40874 6.31561 7.23504

