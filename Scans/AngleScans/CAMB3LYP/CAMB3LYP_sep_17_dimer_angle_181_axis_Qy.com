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
Mg 8.92708 0.80654 0.67116
C 9.47019 -1.91710 -1.60170
C 9.73953 2.90892 -1.88901
C 8.89652 3.19754 2.75673
C 9.27834 -1.48051 3.27726
N 9.51073 0.42538 -1.42796
C 9.62353 -0.64743 -2.19654
C 9.83782 -0.33983 -3.70287
C 10.46767 1.13054 -3.56081
C 9.84358 1.53598 -2.22709
C 11.97867 1.12493 -3.39902
C 8.48953 -0.31212 -4.50497
C 8.64641 -0.32473 -6.02714
H 7.65732 -1.04657 -6.82283
N 8.78687 2.76439 0.37836
C 9.14109 3.44891 -0.72944
C 8.86046 4.86280 -0.43563
C 8.46624 4.97533 0.99498
C 8.72706 3.62352 1.44739
C 8.90018 5.93541 -1.48509
C 7.95312 6.11473 1.82242
O 7.90147 6.04838 3.10139
C 7.60064 7.34449 1.09788
N 9.27630 0.89609 2.71869
C 9.14625 2.05252 3.39225
C 9.17510 1.82100 4.88422
C 9.10863 0.23447 5.07416
C 9.14442 -0.19952 3.59235
C 10.39750 2.37052 5.61917
C 7.97415 -0.31905 5.97928
C 6.63159 -0.08649 5.32442
N 9.12059 -1.33851 0.84260
C 9.28144 -2.08838 2.02831
C 9.48534 -3.49599 1.62894
C 9.45168 -3.46359 0.22669
C 9.26180 -2.12698 -0.21704
C 9.71656 -4.58664 2.51389
C 9.54422 -4.25465 -0.97928
O 9.61702 -5.45791 -1.22544
C 9.68679 -3.27068 -2.21938
C 8.72055 -3.71889 -3.28145
O 7.49387 -3.65958 -3.08562
O 9.39371 -3.93710 -4.42584
C 8.47651 -4.37299 -5.51237
H 10.19621 3.64292 -2.55352
H 8.78395 3.96912 3.52190
H 9.29750 -2.26865 4.03300
H 10.53722 -0.99694 -4.20567
H 10.12458 1.78584 -4.35720
H 12.36820 0.12396 -3.28740
H 12.18175 1.71101 -2.51227
H 12.44861 1.66086 -4.22864
H 8.01457 0.66858 -4.37316
H 7.84880 -1.11694 -4.20106
H 9.66245 -0.60065 -6.33452
H 8.52228 0.69587 -6.39725
H 9.38150 5.58067 -2.40481
H 9.56791 6.72072 -1.11457
H 7.91757 6.37013 -1.68046
H 8.47200 7.70821 0.55853
H 7.17042 8.01195 1.84707
H 6.84293 7.10171 0.34452
H 8.31212 2.26856 5.36177
H 10.02456 -0.15897 5.53904
H 10.13793 3.07409 6.43261
H 11.04475 2.84798 4.89554
H 11.05255 1.61716 6.05341
H 7.95261 0.21810 6.93068
H 8.07983 -1.37437 6.20760
H 6.39799 0.94616 5.60037
H 5.88880 -0.77775 5.71745
H 6.71458 -0.12506 4.23205
H 9.10885 -5.40886 2.16012
H 9.47227 -4.32389 3.54105
H 10.79106 -4.72924 2.37111
H 10.72171 -3.32073 -2.56000
H 8.86073 -5.32472 -5.89221
H 8.43324 -3.52720 -6.22407
H 7.46588 -4.64862 -5.23470

