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
Mg 15.75770 1.41891 1.20386
C 16.53547 3.94055 3.63604
C 15.63835 -0.81555 3.77826
C 15.51777 -0.89866 -0.94990
C 17.04362 3.56194 -1.22503
N 16.06434 1.65310 3.38180
C 16.29424 2.66992 4.19884
C 16.18450 2.32047 5.70728
C 16.47412 0.74388 5.61057
C 15.99874 0.49531 4.18057
C 17.95163 0.39745 5.68976
C 14.75403 2.60735 6.28476
C 14.66020 2.58268 7.81196
H 13.74769 3.51482 8.46839
N 15.12572 -0.45251 1.39606
C 15.12837 -1.20100 2.51902
C 14.58194 -2.51057 2.13094
C 14.41019 -2.52780 0.65279
C 15.04458 -1.27380 0.29907
C 14.20311 -3.56320 3.13187
C 13.79041 -3.51615 -0.28809
O 13.96385 -3.43924 -1.55570
C 13.05190 -4.63017 0.32433
N 16.40484 1.25102 -0.76473
C 16.12374 0.15687 -1.49393
C 16.44698 0.37569 -2.95270
C 16.77860 1.93412 -3.09019
C 16.67195 2.34745 -1.60588
C 17.61311 -0.44323 -3.50579
C 15.96418 2.73689 -4.14136
C 14.51606 2.82326 -3.71607
N 16.46411 3.45994 1.14696
C 16.98338 4.15204 0.03097
C 17.43761 5.47334 0.51084
C 17.17007 5.44930 1.88792
C 16.60865 4.19355 2.24486
C 18.05385 6.48038 -0.28429
C 17.24473 6.19673 3.12261
O 17.55108 7.34988 3.42283
C 16.95400 5.20629 4.33128
C 15.95727 5.86694 5.24373
O 14.79834 6.09502 4.85453
O 16.46746 5.92207 6.48753
C 15.51098 6.55931 7.43140
H 15.79997 -1.63591 4.47810
H 15.35675 -1.62260 -1.75206
H 17.36587 4.32414 -1.93770
H 16.92495 2.79648 6.33889
H 15.86495 0.18631 6.31725
H 18.57355 1.28020 5.67926
H 18.15588 -0.21957 4.82453
H 18.14467 -0.23342 6.56216
H 14.09433 1.76428 6.04234
H 14.37349 3.53935 5.91436
H 15.64865 2.61428 8.28622
H 14.24645 1.61896 8.11868
H 14.59709 -3.33057 4.12890
H 14.72360 -4.48234 2.84181
H 13.12858 -3.75716 3.15292
H 13.71687 -5.18698 0.98026
H 12.60740 -5.17887 -0.50827
H 12.25822 -4.21779 0.95727
H 15.59361 0.14156 -3.57724
H 17.82344 2.10353 -3.38927
H 17.33450 -1.06678 -4.37621
H 18.00694 -1.05850 -2.70774
H 18.48533 0.13696 -3.80212
H 15.97464 2.21978 -5.10398
H 16.34517 3.73863 -4.30976
H 14.09945 1.87349 -4.06458
H 14.02596 3.66861 -4.19485
H 14.42710 2.84115 -2.62363
H 17.60215 7.42143 0.00012
H 17.92591 6.28202 -1.34640
H 19.09450 6.36875 0.03170
H 17.90327 5.01383 4.83271
H 16.03666 7.39527 7.90311
H 15.15951 5.74667 8.09459
H 14.64963 7.06282 7.00831

