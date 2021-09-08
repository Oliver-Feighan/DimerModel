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
Mg 7.87829 0.71387 0.60728
C 5.55766 -1.90219 1.41352
C 5.18738 2.76770 0.18853
C 9.85902 2.91056 -0.54105
C 10.24267 -1.76764 -0.02321
N 5.70354 0.36217 0.80488
C 4.94175 -0.65562 1.17651
C 3.44481 -0.29826 1.37909
C 3.37523 0.99387 0.42808
C 4.83722 1.42976 0.50060
C 3.04752 0.66002 -1.01777
C 3.11368 0.06648 2.86883
C 1.62081 0.12467 3.20009
H 1.18327 -0.31751 4.52126
N 7.64043 2.67276 0.39577
C 6.47578 3.34127 0.26182
C 6.84009 4.75926 0.11770
C 8.32148 4.85459 0.01226
C 8.67070 3.45060 -0.07115
C 5.82990 5.86837 0.17066
C 9.26636 6.01743 -0.02233
O 10.49537 5.87664 -0.35778
C 8.68865 7.34012 0.25740
N 9.70783 0.58871 -0.37222
C 10.38483 1.69815 -0.71682
C 11.79055 1.36489 -1.15648
C 11.99576 -0.18144 -0.80388
C 10.58077 -0.51154 -0.28021
C 12.09666 1.59413 -2.63639
C 13.21639 -0.54452 0.08547
C 13.02249 0.00881 1.47889
N 7.98460 -1.43055 0.84861
C 9.05927 -2.27629 0.49663
C 8.61931 -3.66485 0.74311
C 7.30068 -3.53089 1.20327
C 6.93713 -2.15721 1.22176
C 9.38733 -4.83730 0.49486
C 6.13000 -4.23991 1.66750
O 5.87650 -5.41204 1.94241
C 4.90656 -3.22589 1.70424
C 4.20854 -3.38790 3.02675
O 4.78460 -3.08570 4.08667
O 2.90994 -3.66436 2.80864
C 2.17308 -3.82299 4.09070
H 4.41020 3.43283 -0.18895
H 10.61825 3.63474 -0.84534
H 10.95497 -2.59172 -0.10202
H 2.74715 -1.05190 1.03404
H 2.72788 1.75919 0.84830
H 3.03183 -0.40566 -1.19160
H 3.82199 1.12872 -1.61071
H 2.11010 1.14040 -1.31242
H 3.38738 1.11380 3.05038
H 3.60780 -0.60454 3.54425
H 1.00662 -0.33697 2.41740
H 1.30691 1.17132 3.20466
H 4.80569 5.48379 0.09014
H 5.96656 6.46887 -0.73521
H 5.95646 6.51204 1.04378
H 7.89902 7.54849 -0.46060
H 9.53415 8.03056 0.27470
H 8.21637 7.31368 1.24576
H 12.51689 1.95022 -0.60610
H 12.14582 -0.78915 -1.70831
H 12.94864 2.27959 -2.80479
H 11.20341 1.97303 -3.11500
H 12.30163 0.69226 -3.21067
H 14.12377 -0.08075 -0.30915
H 13.40183 -1.61151 0.14973
H 13.35585 1.04630 1.38102
H 13.63330 -0.53607 2.19586
H 11.96064 0.02800 1.75049
H 9.24754 -5.48753 1.34818
H 10.43815 -4.59952 0.34316
H 8.91023 -5.19238 -0.42251
H 4.25429 -3.46894 0.86457
H 1.69308 -4.80613 4.06379
H 1.51007 -2.94062 4.16537
H 2.75878 -3.89841 4.99940

