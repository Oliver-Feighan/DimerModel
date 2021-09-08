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
Mg 13.14490 1.18774 1.00499
C 12.92406 -1.03498 3.81387
C 11.88821 3.62931 3.02837
C 12.90305 3.04658 -1.55502
C 13.28486 -1.63191 -1.03846
N 12.52565 1.17644 3.12835
C 12.54464 0.29026 4.11256
C 12.21099 0.87513 5.51119
C 11.36633 2.15724 5.04043
C 11.98403 2.36305 3.65885
C 9.88467 1.86983 4.86318
C 13.49719 1.29482 6.30590
C 13.26955 1.58418 7.79123
H 14.31916 1.22205 8.73973
N 12.96210 3.15868 0.86164
C 12.45097 3.99791 1.78669
C 12.52074 5.34200 1.19279
C 12.96234 5.20735 -0.22198
C 12.93908 3.76588 -0.36933
C 12.26192 6.59380 1.97983
C 13.33070 6.21051 -1.27278
O 13.45495 5.88017 -2.50496
C 13.45011 7.60795 -0.83188
N 12.88709 0.77673 -1.01663
C 12.86734 1.76667 -1.92633
C 12.94853 1.21836 -3.33102
C 13.27198 -0.33914 -3.16623
C 13.23179 -0.44435 -1.62587
C 11.69278 1.38946 -4.18554
C 14.52229 -0.88100 -3.91174
C 15.77799 -0.29666 -3.30581
N 13.29849 -0.94607 1.30718
C 13.31561 -1.94924 0.31348
C 13.31534 -3.25204 1.01000
C 13.27451 -2.91394 2.37107
C 13.23060 -1.50057 2.51233
C 13.30164 -4.53019 0.38369
C 13.24784 -3.43152 3.72018
O 13.35242 -4.54905 4.22409
C 12.89206 -2.24091 4.71127
C 13.86317 -2.28608 5.85906
O 15.07354 -2.07035 5.67178
O 13.17705 -2.36035 7.01409
C 14.09665 -2.39642 8.18236
H 11.29006 4.40413 3.50874
H 12.93085 3.64343 -2.46961
H 13.42660 -2.55632 -1.60242
H 11.59929 0.23608 6.13666
H 11.56294 3.01549 5.67770
H 9.66273 0.81834 4.96926
H 9.63615 2.21029 3.86650
H 9.29635 2.48726 5.54801
H 13.81855 2.28845 5.96801
H 14.27025 0.55944 6.19535
H 12.29541 1.21834 8.13787
H 13.21391 2.66670 7.92933
H 11.79747 6.37118 2.94840
H 11.49845 7.16102 1.43647
H 13.15348 7.21513 2.08817
H 12.50707 7.93150 -0.39777
H 13.80678 8.15996 -1.70363
H 14.19852 7.66010 -0.03327
H 13.75317 1.68784 -3.88351
H 12.45298 -0.96774 -3.54562
H 11.87862 1.93451 -5.13024
H 10.94393 1.89880 -3.59365
H 11.18600 0.46359 -4.45226
H 14.50618 -0.56479 -4.95764
H 14.59460 -1.96345 -3.90364
H 15.86022 0.67651 -3.79903
H 16.63832 -0.92515 -3.52688
H 15.64845 -0.11258 -2.23297
H 14.01251 -5.14636 0.91801
H 13.55197 -4.45801 -0.67267
H 12.25704 -4.81261 0.53963
H 11.86228 -2.38523 5.04036
H 13.84806 -3.29399 8.75704
H 13.97183 -1.42202 8.69100
H 15.15040 -2.55612 7.98600

