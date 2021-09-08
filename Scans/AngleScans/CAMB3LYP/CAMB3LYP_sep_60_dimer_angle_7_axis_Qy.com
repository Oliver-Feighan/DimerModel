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
Mg 31.53977 2.84368 2.40136
C 31.64370 0.65080 5.24029
C 30.57939 5.32210 4.53905
C 31.00326 4.68180 -0.12864
C 31.38482 0.00329 0.38794
N 31.19176 2.86127 4.58559
C 31.32236 1.98415 5.56944
C 31.17463 2.58688 6.99228
C 30.29464 3.87668 6.61688
C 30.73677 4.06054 5.16658
C 28.79877 3.60920 6.62994
C 32.55583 2.99524 7.61480
C 32.52023 3.30177 9.11364
H 33.67567 2.93330 9.92707
N 31.36657 4.81565 2.26017
C 30.98673 5.67084 3.23261
C 30.99924 7.00816 2.61975
C 31.25799 6.85386 1.16237
C 31.19732 5.41157 1.03512
C 30.85784 8.27090 3.41907
C 31.50485 7.84168 0.06261
O 31.46911 7.49805 -1.17168
C 31.69716 9.24130 0.46951
N 31.02490 2.41756 0.43269
C 30.90428 3.39913 -0.47830
C 30.80129 2.83657 -1.87591
C 31.12219 1.27614 -1.73567
C 31.27421 1.18595 -0.20134
C 29.45062 3.01791 -2.56814
C 32.26177 0.70920 -2.62594
C 33.59121 1.28086 -2.18871
N 31.70177 0.71079 2.70561
C 31.58076 -0.30178 1.72880
C 31.65062 -1.59786 2.43430
C 31.78539 -1.24648 3.78591
C 31.77829 0.16864 3.91585
C 31.54150 -2.88149 1.82888
C 31.92136 -1.75097 5.13339
O 32.07352 -2.86514 5.63260
C 31.70857 -0.54609 6.14797
C 32.81535 -0.59463 7.16542
O 33.99540 -0.39829 6.82552
O 32.27866 -0.64810 8.39809
C 33.33703 -0.68661 9.44221
H 30.05657 6.11000 5.08196
H 30.92399 5.26960 -1.04607
H 31.44241 -0.92832 -0.17903
H 30.63783 1.96268 7.69654
H 30.58102 4.73790 7.21491
H 28.57798 2.56209 6.77468
H 28.43169 3.94390 5.66855
H 28.30927 4.24151 7.37621
H 32.84540 3.98088 7.22827
H 33.29909 2.24769 7.41640
H 31.59252 2.95339 9.58368
H 32.49671 4.38623 9.24559
H 30.51569 8.06414 4.44063
H 30.03981 8.84406 2.96942
H 31.76410 8.88019 3.40790
H 30.82042 9.58258 1.01478
H 31.94891 9.77988 -0.44614
H 32.54049 9.29003 1.16739
H 31.53639 3.28911 -2.53008
H 30.25379 0.65600 -2.00241
H 29.52365 3.55134 -3.53468
H 28.78878 3.54361 -1.89276
H 28.90215 2.09706 -2.75894
H 32.11874 1.01581 -3.66502
H 32.32017 -0.37406 -2.61495
H 33.62380 2.24807 -2.69913
H 34.40860 0.63788 -2.50887
H 33.59977 1.47684 -1.11021
H 32.30557 -3.50290 2.27670
H 31.65823 -2.82284 0.74874
H 30.52107 -3.14721 2.11765
H 30.72640 -0.67233 6.60512
H 33.15064 -1.57505 10.05346
H 33.28995 0.29422 9.95163
H 34.35560 -0.86345 9.11711

