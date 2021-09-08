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
Mg 31.54063 2.83499 2.38962
C 31.36529 4.93070 -0.51835
C 31.05484 0.17108 0.31545
C 31.24020 1.00941 4.96731
C 30.87097 5.67016 4.30248
N 31.26047 2.69640 0.19997
C 31.27396 3.55365 -0.80963
C 31.26538 2.89318 -2.21436
C 30.59654 1.48870 -1.81564
C 31.02249 1.42234 -0.35043
C 29.07803 1.50796 -1.87120
C 32.71121 2.69752 -2.79170
C 32.76721 2.34520 -4.27982
H 33.86946 2.87314 -5.07907
N 31.68617 0.86623 2.59309
C 31.47723 -0.06781 1.64167
C 31.68985 -1.36648 2.29925
C 31.87982 -1.12909 3.75605
C 31.58223 0.28709 3.83337
C 31.77751 -2.65839 1.53985
C 32.25341 -2.03049 4.89364
O 32.12831 -1.66103 6.11458
C 32.68162 -3.39134 4.53841
N 30.90929 3.22927 4.33083
C 30.92461 2.26849 5.27136
C 30.69312 2.84768 6.64662
C 30.76007 4.43474 6.46124
C 30.93769 4.50322 4.92855
C 29.37126 2.46885 7.31412
C 31.76746 5.20608 7.35727
C 33.18364 4.84640 6.96939
N 31.36241 4.95581 2.01758
C 31.05166 5.96343 2.95685
C 30.92953 7.23214 2.20994
C 31.15680 6.86774 0.87426
C 31.38321 5.46723 0.79182
C 30.59673 8.49804 2.76928
C 31.24615 7.34767 -0.48603
O 31.22905 8.45651 -1.01886
C 31.25993 6.09481 -1.46406
C 32.37176 6.29329 -2.45751
O 33.55821 6.30229 -2.08486
O 31.86772 6.22206 -3.70300
C 32.93417 6.40195 -4.72385
H 30.68210 -0.70726 -0.21226
H 31.23217 0.44389 5.90197
H 30.76083 6.61500 4.83898
H 30.65397 3.40035 -2.95105
H 31.03537 0.66850 -2.37782
H 28.69418 2.50027 -2.05595
H 28.74380 1.14636 -0.90757
H 28.71849 0.78251 -2.60657
H 33.14623 1.78419 -2.36588
H 33.31743 3.56188 -2.60193
H 31.80856 2.52348 -4.78197
H 32.92378 1.26803 -4.37574
H 31.43462 -2.54047 0.50451
H 31.05134 -3.34385 1.99006
H 32.77004 -3.11100 1.59179
H 31.88732 -3.88707 3.98542
H 32.99214 -3.85436 5.47706
H 33.54059 -3.32236 3.86163
H 31.47364 2.54069 7.33193
H 29.79547 4.91244 6.68747
H 29.50328 1.98325 8.29944
H 28.82247 1.82249 6.64212
H 28.67547 3.29307 7.46153
H 31.64752 4.91103 8.40261
H 31.64944 6.28355 7.31119
H 33.35879 3.91297 7.51263
H 33.87660 5.62332 7.28600
H 33.25366 4.62271 5.89855
H 31.26179 9.22234 2.31807
H 30.69161 8.49114 3.85313
H 29.55505 8.58503 2.44896
H 30.28328 6.04552 -1.94708
H 32.62289 7.22964 -5.36863
H 33.06108 5.41199 -5.20095
H 33.90114 6.75217 -4.38214

