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
Mg 8.94159 0.80883 0.67252
C 8.36874 4.28364 1.36291
C 8.24467 1.49680 -2.59483
C 8.99662 -2.31276 0.10679
C 8.46924 0.19863 4.07109
N 8.39200 2.68631 -0.35960
C 8.26472 3.95613 -0.00506
C 8.09147 4.94111 -1.19223
C 7.50135 3.91187 -2.27449
C 8.10834 2.61746 -1.73715
C 5.98771 3.78690 -2.22329
C 9.45676 5.55248 -1.66577
C 9.33200 6.74970 -2.61068
H 10.31659 7.82322 -2.50772
N 9.14433 -0.18368 -1.03405
C 8.83362 0.25450 -2.27208
C 9.14853 -0.85688 -3.18303
C 9.51530 -2.04457 -2.36475
C 9.20609 -1.55507 -1.03626
C 9.16206 -0.70600 -4.67651
C 10.04331 -3.40141 -2.72027
O 10.06581 -4.35944 -1.86910
C 10.44640 -3.60197 -4.11984
N 8.55161 -0.84764 1.86682
C 8.70055 -2.09536 1.38835
C 8.63332 -3.11692 2.49844
C 8.65001 -2.27541 3.85828
C 8.63327 -0.84666 3.27194
C 7.41189 -4.03577 2.48383
C 9.74856 -2.63588 4.89554
C 11.11089 -2.26740 4.35362
N 8.68278 2.01291 2.44770
C 8.47521 1.55393 3.76692
C 8.23939 2.74398 4.60989
C 8.30380 3.81966 3.71142
C 8.54149 3.33276 2.39767
C 7.95808 2.73256 6.00519
C 8.21419 5.25554 3.57391
O 8.11202 6.19855 4.35752
C 8.12643 5.61057 2.02710
C 9.10174 6.72152 1.74945
O 10.32517 6.53312 1.86913
O 8.44707 7.75109 1.18220
C 9.37424 8.87340 0.87792
H 7.82356 1.55754 -3.59876
H 9.11506 -3.39198 -0.01464
H 8.41132 0.09972 5.15715
H 7.38429 5.74309 -1.01759
H 7.88004 4.12586 -3.27056
H 5.56723 4.32366 -1.38590
H 5.78278 2.72798 -2.13533
H 5.55116 4.10024 -3.17596
H 9.95685 4.83824 -2.33248
H 10.06743 5.81797 -0.82481
H 8.31520 7.16070 -2.62322
H 9.49331 6.40234 -3.63408
H 8.69040 0.23344 -4.99019
H 8.50949 -1.48642 -5.08266
H 10.16071 -0.82322 -5.10268
H 9.59754 -3.40857 -4.77139
H 10.87963 -4.60331 -4.15855
H 11.21266 -2.86073 -4.37288
H 9.49847 -3.76820 2.47757
H 7.71345 -2.39448 4.42264
H 7.67426 -5.11008 2.45087
H 6.79420 -3.76902 1.63657
H 6.72633 -3.90597 3.31951
H 9.76532 -3.71450 5.07006
H 9.60770 -2.15603 5.85829
H 11.36821 -3.12703 3.72759
H 11.82494 -2.14228 5.16511
H 11.05009 -1.38977 3.69962
H 8.54925 3.52283 6.44846
H 8.18793 1.76623 6.44921
H 6.88320 2.93166 5.98907
H 7.09798 5.90660 1.81677
H 8.97094 9.76352 1.37059
H 9.45699 8.90195 -0.22484
H 10.37038 8.82972 1.30249

