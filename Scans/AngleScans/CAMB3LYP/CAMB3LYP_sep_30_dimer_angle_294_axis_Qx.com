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
Mg 15.77239 1.41558 1.20543
C 13.30179 0.62280 3.68468
C 17.16023 3.54130 3.48356
C 17.68837 2.46974 -1.09353
C 13.57948 0.14393 -1.18617
N 15.21275 1.92886 3.28284
C 14.29116 1.51849 4.14110
C 14.48524 2.03852 5.59074
C 15.34569 3.35378 5.26121
C 15.97804 2.90855 3.94417
C 14.49033 4.58018 4.98983
C 15.29307 1.03239 6.48356
C 15.25459 1.33014 7.98410
H 15.24056 0.20639 8.91655
N 17.48403 2.41378 1.31716
C 17.89653 3.21153 2.32449
C 19.21304 3.72611 1.91686
C 19.47050 3.29798 0.51503
C 18.18723 2.71545 0.17742
C 20.12048 4.47150 2.85194
C 20.66866 3.40766 -0.37860
O 20.59245 3.18123 -1.63786
C 21.91271 3.87843 0.24768
N 15.57120 1.51790 -0.86130
C 16.57542 1.97109 -1.63206
C 16.30906 1.69552 -3.09272
C 15.02454 0.74333 -3.12301
C 14.70311 0.70390 -1.61287
C 16.05236 2.92390 -3.96540
C 15.16453 -0.60332 -3.88433
C 16.13309 -1.50862 -3.15788
N 13.88067 0.37208 1.22834
C 13.12316 -0.04525 0.11205
C 11.85885 -0.60758 0.62949
C 11.95944 -0.45619 2.02069
C 13.19014 0.17556 2.34593
C 10.79151 -1.12021 -0.16052
C 11.29970 -0.67750 3.28736
O 10.27193 -1.25300 3.64264
C 12.08370 0.12589 4.41266
C 12.29399 -0.79737 5.58138
O 13.02094 -1.80042 5.47165
O 11.78190 -0.23009 6.68877
C 11.97595 -1.10322 7.87696
H 17.53908 4.39998 4.03855
H 18.36044 2.71448 -1.91930
H 12.88631 -0.36208 -1.86154
H 13.56977 2.31150 6.10180
H 16.11145 3.52360 6.01364
H 13.43911 4.33675 4.94708
H 14.82758 4.96929 4.03804
H 14.70417 5.35952 5.72687
H 16.36570 1.16606 6.29293
H 14.97254 0.02409 6.30690
H 14.47030 2.05251 8.24102
H 16.18389 1.83335 8.26197
H 19.59369 4.77915 3.76367
H 20.38334 5.41376 2.35907
H 21.03734 3.92221 3.07625
H 21.74820 4.86093 0.68367
H 22.67667 3.80517 -0.52877
H 22.16454 3.20526 1.07475
H 17.14334 1.17472 -3.54641
H 14.16858 1.22932 -3.61390
H 16.74836 3.00413 -4.82169
H 16.11453 3.80524 -3.34098
H 15.04523 2.99932 -4.37187
H 15.58425 -0.43417 -4.87901
H 14.22235 -1.12461 -4.01673
H 17.10822 -1.16352 -3.51439
H 15.96017 -2.54882 -3.42645
H 16.09756 -1.33568 -2.07599
H 10.46388 -2.03518 0.31488
H 11.10543 -1.29519 -1.18746
H 10.07853 -0.29492 -0.08419
H 11.47725 0.99230 4.67934
H 10.98775 -1.26731 8.31767
H 12.72864 -0.59038 8.50482
H 12.30151 -2.12204 7.70245

