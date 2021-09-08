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
Mg 13.14009 1.19130 1.00424
C 14.62832 -0.33664 3.89035
C 10.63969 2.27813 3.05398
C 11.57679 2.16033 -1.58120
C 15.08883 -0.93361 -0.95353
N 12.79089 0.93988 3.17383
C 13.45600 0.39227 4.17985
C 12.87315 0.70751 5.58359
C 11.35935 1.02277 5.15015
C 11.60615 1.47702 3.71299
C 10.47376 -0.21169 5.11548
C 13.55294 1.95682 6.24626
C 13.25370 2.13706 7.73613
H 14.30533 2.67011 8.59776
N 11.64785 2.48600 0.81571
C 10.74094 2.82494 1.75588
C 9.84280 3.79862 1.11599
C 10.19414 3.88266 -0.32769
C 11.16063 2.80741 -0.42680
C 8.82969 4.59768 1.88315
C 9.72690 4.77385 -1.43839
O 9.99026 4.51468 -2.66571
C 8.87319 5.90749 -1.05428
N 13.14723 0.54628 -0.97218
C 12.41356 1.17502 -1.90724
C 12.78784 0.71361 -3.29555
C 14.09919 -0.18305 -3.11132
C 14.20931 -0.15687 -1.57101
C 11.72188 -0.09088 -4.03913
C 15.34614 0.21340 -3.94835
C 15.88252 1.54555 -3.47607
N 14.72988 -0.22821 1.35865
C 15.38788 -1.02901 0.39956
C 16.31257 -1.91663 1.13408
C 16.11002 -1.58404 2.48198
C 15.11382 -0.57554 2.58193
C 17.15301 -2.90716 0.55216
C 16.50478 -1.86428 3.84372
O 17.36996 -2.56208 4.37125
C 15.47221 -1.15859 4.82462
C 16.25790 -0.43206 5.88167
O 16.97996 0.53450 5.58001
O 15.86126 -0.85656 7.09528
C 16.60422 -0.15683 8.17704
H 9.69450 2.47338 3.56130
H 11.14733 2.53559 -2.51298
H 15.80188 -1.55593 -1.49846
H 12.89528 -0.12113 6.28119
H 10.94044 1.83403 5.73984
H 11.03933 -1.11783 5.27414
H 10.01626 -0.21807 4.13482
H 9.65268 -0.10678 5.83042
H 13.08916 2.86920 5.84975
H 14.61407 1.94074 6.09017
H 12.81308 1.23607 8.18000
H 12.47603 2.89700 7.84443
H 8.68766 4.20090 2.89597
H 7.86252 4.44297 1.39273
H 9.05476 5.66629 1.89154
H 7.98557 5.53581 -0.54788
H 8.71504 6.47815 -1.97146
H 9.41525 6.52332 -0.32796
H 13.02534 1.55686 -3.93243
H 13.91997 -1.23033 -3.39581
H 11.44136 0.35203 -5.01341
H 10.85447 -0.18184 -3.39892
H 11.97826 -1.13157 -4.23004
H 15.07180 0.34379 -4.99806
H 16.14223 -0.52278 -3.91261
H 15.25249 2.26684 -4.00514
H 16.92877 1.65731 -3.75337
H 15.70881 1.68143 -2.40229
H 18.11522 -2.82465 1.03984
H 17.23913 -2.77313 -0.52413
H 16.59567 -3.81172 0.80987
H 14.83833 -1.93812 5.24901
H 17.06515 -0.92981 8.79971
H 15.86671 0.50867 8.66381
H 17.47002 0.42945 7.89222

