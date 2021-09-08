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
Mg 13.13098 1.18506 0.99099
C 13.51162 -1.55351 -1.29690
C 13.79520 3.26901 -1.62637
C 13.24782 3.58642 3.06155
C 13.62977 -1.09162 3.58207
N 13.57893 0.78945 -1.13858
C 13.63591 -0.28789 -1.90698
C 13.75719 0.01078 -3.42542
C 14.40471 1.47730 -3.33107
C 13.86845 1.89371 -1.96297
C 15.92282 1.46172 -3.26448
C 12.36137 0.04415 -4.14138
C 12.42219 0.02295 -5.67031
H 11.38016 -0.69571 -6.39843
N 12.98596 3.14240 0.69711
C 13.27451 3.81893 -0.43440
C 13.02253 5.23621 -0.13110
C 12.71978 5.35856 1.32083
C 12.99931 4.00717 1.76318
C 13.00350 6.30334 -1.18670
C 12.26743 6.50565 2.17276
O 12.29581 6.44594 3.45278
C 11.87849 7.73431 1.46522
N 13.60880 1.28217 3.01203
C 13.52921 2.44279 3.68624
C 13.65019 2.21840 5.17466
C 13.58500 0.63333 5.37688
C 13.52463 0.19183 3.89809
C 14.92008 2.76279 5.82843
C 12.50591 0.09236 6.35443
C 11.12646 0.33126 5.78396
N 13.32028 -0.96045 1.16141
C 13.55023 -1.70562 2.33866
C 13.71905 -3.11658 1.93480
C 13.59754 -3.09083 0.53728
C 13.38924 -1.75509 0.09921
C 13.99801 -4.20450 2.80930
C 13.60872 -3.88844 -0.66787
O 13.65772 -5.09338 -0.91168
C 13.67977 -2.91161 -1.91973
C 12.64566 -3.35813 -2.91659
O 11.43414 -3.28912 -2.64442
O 13.24406 -3.58675 -4.09982
C 12.25744 -4.02142 -5.12425
H 14.21420 3.99647 -2.32216
H 13.18882 4.36253 3.82814
H 13.69102 -1.87615 4.33932
H 14.41912 -0.65376 -3.96764
H 14.01671 2.13111 -4.10783
H 16.31178 0.45856 -3.17220
H 16.18521 2.05068 -2.39539
H 16.34333 1.99021 -4.12483
H 11.90232 1.02885 -3.98525
H 11.73555 -0.75458 -3.79353
H 13.41500 -0.26171 -6.03943
H 12.28199 1.04258 -6.03734
H 13.42364 5.94068 -2.13293
H 13.69853 7.08569 -0.86306
H 12.01354 6.74409 -1.32229
H 12.71668 8.08916 0.87027
H 11.50075 8.40849 2.23638
H 11.07329 7.49324 0.76226
H 12.82200 2.67444 5.70308
H 14.52564 0.23566 5.78541
H 14.71693 3.47218 6.65279
H 15.52380 3.23207 5.06303
H 15.59598 2.00692 6.22469
H 12.54786 0.63431 7.30241
H 12.61854 -0.96255 6.58129
H 10.91770 1.36689 6.06851
H 10.40515 -0.35275 6.22656
H 11.14036 0.28674 4.68875
H 13.36368 -5.02409 2.49880
H 13.82055 -3.93498 3.84835
H 15.06041 -4.35545 2.60008
H 14.69086 -2.97071 -2.32440
H 12.61054 -4.97771 -5.52238
H 12.17528 -3.17884 -5.83634
H 11.26441 -4.28847 -4.78217

