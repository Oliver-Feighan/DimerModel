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
Mg 7.89258 0.71427 0.60039
C 9.25912 -1.40744 3.15185
C 8.19185 3.30254 2.80173
C 6.28144 2.51741 -1.45396
C 6.66207 -2.16113 -0.93702
N 8.64343 0.80978 2.67869
C 9.19467 -0.04839 3.52369
C 9.77846 0.59870 4.80822
C 8.88335 1.93174 4.83408
C 8.57786 2.05000 3.34232
C 7.56968 1.76020 5.57852
C 11.30477 0.93434 4.66778
C 12.01050 1.28075 5.98068
H 13.39789 0.86060 6.15704
N 7.75827 2.68899 0.45545
C 7.93252 3.59117 1.44387
C 7.70547 4.90888 0.83047
C 7.22153 4.70123 -0.56144
C 7.04460 3.26292 -0.56719
C 8.02234 6.19809 1.53115
C 6.94927 5.64292 -1.69515
O 6.30725 5.27079 -2.74021
C 7.37450 7.03746 -1.50559
N 6.47336 0.27147 -0.85292
C 5.97054 1.23508 -1.64446
C 5.18137 0.64463 -2.78856
C 5.46222 -0.92893 -2.73725
C 6.33202 -0.98930 -1.46242
C 3.67363 0.89359 -2.75510
C 6.00464 -1.58967 -4.03412
C 7.40330 -1.09257 -4.32054
N 8.08883 -1.41613 0.90198
C 7.46759 -2.44363 0.15876
C 7.81319 -3.72306 0.81148
C 8.59880 -3.34648 1.91125
C 8.71658 -1.93090 1.95336
C 7.36990 -5.01223 0.40203
C 9.34642 -3.82394 3.05221
O 9.67230 -4.93238 3.47517
C 9.70227 -2.58242 3.97870
C 11.15944 -2.67509 4.34008
O 12.03600 -2.56283 3.46507
O 11.28280 -2.66282 5.67985
C 12.71094 -2.74172 6.08731
H 8.03076 4.13568 3.48643
H 5.79463 3.08539 -2.25013
H 6.39840 -3.10872 -1.41160
H 9.62192 0.02806 5.71577
H 9.45983 2.78807 5.17429
H 7.40111 0.73327 5.86684
H 6.79895 2.09282 4.89559
H 7.52918 2.44119 6.43347
H 11.41411 1.88940 4.13802
H 11.82675 0.13628 4.17659
H 11.41091 1.00413 6.85642
H 12.10058 2.36754 6.04998
H 8.20730 6.03967 2.60077
H 7.11455 6.81033 1.50024
H 8.83600 6.74812 1.05357
H 6.88565 7.44756 -0.62507
H 7.17595 7.53537 -2.45666
H 8.45116 7.05039 -1.30252
H 5.52812 1.03270 -3.73844
H 4.54704 -1.49931 -2.52057
H 3.29398 1.39640 -3.66456
H 3.44355 1.47735 -1.87372
H 3.06193 0.00483 -2.60972
H 5.39116 -1.30128 -4.89123
H 6.01414 -2.67382 -3.99511
H 7.22725 -0.14275 -4.83429
H 7.93585 -1.79413 -4.95949
H 7.93994 -0.87000 -3.39095
H 8.22751 -5.66924 0.45983
H 6.95308 -4.98877 -0.60270
H 6.60527 -5.20516 1.15925
H 9.05843 -2.63431 4.85754
H 12.80454 -3.60060 6.75905
H 12.95814 -1.74717 6.50384
H 13.43725 -2.99113 5.32264

