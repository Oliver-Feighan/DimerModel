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
Mg 15.76888 1.41266 1.20542
C 13.10158 1.99749 3.53406
C 17.92612 2.38922 3.66046
C 18.09700 1.34075 -0.94914
C 13.39773 1.59199 -1.34236
N 15.43622 2.08319 3.28676
C 14.38429 2.20675 4.08207
C 14.73119 2.49684 5.56702
C 16.18110 3.15273 5.35153
C 16.56170 2.47748 4.03554
C 16.13432 4.65420 5.12191
C 14.81513 1.18774 6.42797
C 14.84391 1.41396 7.94116
H 14.17003 0.44629 8.80253
N 17.73736 1.33051 1.44395
C 18.44573 1.75063 2.51306
C 19.85662 1.48899 2.18891
C 19.93612 1.03221 0.77478
C 18.56559 1.24125 0.35281
C 20.95862 1.60115 3.20204
C 21.06250 0.50767 -0.06326
O 20.95949 0.39590 -1.33594
C 22.32134 0.21681 0.63833
N 15.78922 1.66953 -0.85757
C 16.92795 1.53500 -1.55980
C 16.65165 1.49011 -3.04372
C 15.06224 1.37900 -3.18162
C 14.67167 1.47286 -1.69030
C 17.14952 2.68967 -3.84968
C 14.50959 0.19192 -4.01713
C 14.79372 -1.11356 -3.30976
N 13.61566 1.54935 1.09142
C 12.82728 1.63857 -0.07674
C 11.42764 1.82858 0.35608
C 11.50242 1.86002 1.75676
C 12.85580 1.72109 2.16733
C 10.30596 1.99414 -0.50447
C 10.74549 1.98980 2.98106
O 9.54866 2.04649 3.26055
C 11.76211 2.21142 4.18264
C 11.36899 1.28500 5.30054
O 11.45199 0.05229 5.15859
O 11.16868 2.00482 6.41952
C 10.78743 1.12877 7.55925
H 18.66831 2.89245 4.28079
H 18.84803 1.21074 -1.73178
H 12.58729 1.55856 -2.07362
H 14.07264 3.20347 6.05767
H 16.86781 2.86150 6.14209
H 15.12139 5.01537 5.02362
H 16.68877 4.82963 4.20930
H 16.68316 5.17375 5.91259
H 15.80265 0.72962 6.28767
H 14.01732 0.51564 6.17818
H 14.55285 2.43651 8.21067
H 15.87774 1.33014 8.28496
H 20.61984 2.11597 4.10955
H 21.71629 2.26883 2.77778
H 21.42248 0.63866 3.42815
H 22.67981 1.12003 1.12633
H 22.97644 -0.23207 -0.11082
H 22.11933 -0.51086 1.43224
H 17.10549 0.61652 -3.49534
H 14.63324 2.26347 -3.67509
H 17.83470 2.40909 -4.67180
H 17.63252 3.38010 -3.17104
H 16.36806 3.30683 -4.28970
H 15.01864 0.13895 -4.98265
H 13.44572 0.26302 -4.21745
H 15.82306 -1.33607 -3.60649
H 14.10920 -1.88910 -3.64741
H 14.78559 -0.98145 -2.22158
H 9.50948 1.38479 -0.09851
H 10.54393 1.70901 -1.52713
H 10.14184 3.07054 -0.40606
H 11.69727 3.25947 4.47754
H 9.83841 1.50830 7.95042
H 11.65476 1.13741 8.24585
H 10.52810 0.10046 7.33611

