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
Mg 13.14493 1.18690 1.00549
C 12.97724 -1.37334 3.51470
C 11.89093 3.34181 3.33324
C 12.85617 3.35013 -1.29765
C 13.28543 -1.35201 -1.37776
N 12.55117 0.90249 3.11708
C 12.58972 -0.10011 3.98175
C 12.26761 0.30091 5.44633
C 11.40610 1.62327 5.14965
C 12.00544 2.00740 3.79839
C 9.92503 1.34522 4.95354
C 13.55943 0.63058 6.27369
C 13.34697 0.72863 7.78597
H 14.41097 0.26097 8.67009
N 12.94303 3.15825 1.11314
C 12.43557 3.86928 2.14183
C 12.48637 5.27800 1.72099
C 12.91225 5.32675 0.29586
C 12.89997 3.91506 -0.03136
C 12.22590 6.41824 2.66195
C 13.25918 6.45773 -0.62430
O 13.37164 6.28615 -1.88951
C 13.37148 7.78984 -0.01244
N 12.86668 1.03058 -1.04889
C 12.82734 2.12673 -1.82659
C 12.89662 1.76016 -3.28988
C 13.23577 0.19770 -3.32577
C 13.21487 -0.10066 -1.81054
C 11.62931 2.02439 -4.10259
C 14.48184 -0.23331 -4.14687
C 15.73946 0.28311 -3.48575
N 13.32096 -0.96628 1.03523
C 13.33510 -1.83639 -0.07688
C 13.35465 -3.21633 0.45021
C 13.32706 -3.05239 1.84335
C 13.27235 -1.66850 2.16172
C 13.34477 -4.40571 -0.33171
C 13.32105 -3.73566 3.11686
O 13.44150 -4.90652 3.47505
C 12.96660 -2.68275 4.25358
C 13.95169 -2.86183 5.37610
O 15.15778 -2.61186 5.20448
O 13.28007 -3.08770 6.51990
C 14.21381 -3.26086 7.66444
H 11.29173 4.04393 3.91363
H 12.86779 4.05744 -2.13015
H 13.42860 -2.19672 -2.05500
H 11.66908 -0.41792 5.99297
H 11.60270 2.39661 5.88766
H 9.71367 0.28652 4.92890
H 9.66165 1.80567 4.01032
H 9.33948 1.86563 5.71686
H 13.86795 1.66204 6.06003
H 14.33759 -0.07710 6.06324
H 12.38030 0.31215 8.09426
H 13.28342 1.78459 8.05972
H 11.77502 6.07091 3.59975
H 11.45103 7.04138 2.20246
H 13.11317 7.03015 2.83802
H 12.43086 8.04657 0.46900
H 13.71284 8.45066 -0.81161
H 14.12887 7.74889 0.77832
H 13.69044 2.30358 -3.78752
H 12.41789 -0.38661 -3.77242
H 11.79905 2.68572 -4.97317
H 10.88310 2.44760 -3.44335
H 11.12759 1.13423 -4.47820
H 14.45047 0.21165 -5.14446
H 14.56380 -1.30740 -4.27576
H 15.80720 1.31134 -3.85350
H 16.60261 -0.30376 -3.79332
H 15.62110 0.32956 -2.39695
H 14.06738 -5.07681 0.11322
H 13.58185 -4.19879 -1.37322
H 12.30464 -4.71620 -0.20135
H 11.94213 -2.87784 4.57292
H 13.98003 -4.22604 8.12429
H 14.08646 -2.35944 8.29290
H 15.26652 -3.38379 7.43827

