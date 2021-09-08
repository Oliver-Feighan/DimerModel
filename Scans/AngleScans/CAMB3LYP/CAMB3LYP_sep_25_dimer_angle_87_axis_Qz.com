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
Mg 13.14597 1.17728 1.00321
C 13.35883 -1.74911 -1.06319
C 12.25325 -0.88315 3.57075
C 12.51309 3.76430 2.72779
C 12.92981 3.05802 -1.92284
N 12.87584 -1.01513 1.11559
C 13.04455 -2.03103 0.28262
C 12.94535 -3.43058 0.94669
C 12.04788 -3.03092 2.21700
C 12.43745 -1.55911 2.33818
C 10.55441 -3.10905 1.94724
C 14.34634 -3.98497 1.38476
C 14.36298 -5.46932 1.75648
H 15.54796 -6.25627 1.42635
N 12.96068 1.39811 2.96674
C 12.61257 0.45098 3.86286
C 12.59840 1.12166 5.17212
C 12.80570 2.57923 4.95477
C 12.74578 2.64117 3.50816
C 12.48094 0.37367 6.46832
C 13.00967 3.72916 5.89405
O 12.93128 4.94522 5.49669
C 13.21124 3.39085 7.31052
N 12.56293 3.10584 0.49032
C 12.40640 4.05391 1.43085
C 12.25579 5.42107 0.80749
C 12.58715 5.22441 -0.74464
C 12.79400 3.69402 -0.76733
C 10.88069 6.07188 0.95557
C 13.69639 6.12923 -1.34761
C 15.03846 5.76511 -0.75461
N 13.32645 0.78623 -1.11400
C 13.17445 1.71297 -2.16857
C 13.27409 0.95454 -3.43241
C 13.45559 -0.37478 -3.02196
C 13.44797 -0.44301 -1.60253
C 13.14816 1.49908 -4.74150
C 13.64125 -1.73727 -3.46675
O 13.81513 -2.27888 -4.55770
C 13.46032 -2.70526 -2.21909
C 14.60276 -3.68369 -2.22078
O 15.76925 -3.29363 -2.03704
O 14.11048 -4.93590 -2.22140
C 15.20546 -5.94225 -2.21200
H 11.74722 -1.40949 4.38056
H 12.39909 4.70310 3.27473
H 12.97058 3.58550 -2.87824
H 12.43623 -4.18013 0.35278
H 12.35222 -3.58019 3.10416
H 10.34272 -3.30713 0.90701
H 10.15213 -2.14769 2.23880
H 10.08947 -3.84401 2.61056
H 14.61838 -3.54571 2.35316
H 15.08479 -3.79297 0.63080
H 13.45386 -5.98695 1.42707
H 14.34022 -5.55458 2.84562
H 12.17611 -0.66750 6.30572
H 11.64535 0.81840 7.01955
H 13.38401 0.44377 7.07841
H 12.35322 2.83004 7.67350
H 13.42828 4.33754 7.80906
H 14.07868 2.72632 7.39149
H 12.96550 6.12019 1.23250
H 11.71208 5.43266 -1.37767
H 10.91735 7.06277 1.44635
H 10.24139 5.39688 1.50892
H 10.32913 6.20261 0.02611
H 13.51536 7.17493 -1.08704
H 13.75909 6.07306 -2.42923
H 15.04934 6.31808 0.18942
H 15.84627 6.08589 -1.40927
H 15.08468 4.69717 -0.51163
H 13.92994 1.05214 -5.34115
H 13.22617 2.58423 -4.72992
H 12.13963 1.16274 -4.99647
H 12.49550 -3.20226 -2.32725
H 15.04417 -6.59794 -3.07324
H 15.17297 -6.40973 -1.20993
H 16.21246 -5.58901 -2.40077

