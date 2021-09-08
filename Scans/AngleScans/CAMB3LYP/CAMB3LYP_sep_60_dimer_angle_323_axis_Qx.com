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
Mg 31.54003 2.83882 2.40284
C 29.85218 1.08749 5.04160
C 31.80969 5.48721 4.53542
C 32.62300 4.56150 -0.03170
C 30.14931 0.54435 0.17865
N 30.87803 3.12858 4.49330
C 30.30139 2.37147 5.41451
C 30.27260 2.99660 6.83505
C 30.37772 4.54408 6.41895
C 31.09809 4.38907 5.08115
C 29.02700 5.19125 6.16192
C 31.49790 2.55217 7.70862
C 31.37495 2.87360 9.19970
H 31.93961 1.93288 10.16337
N 32.55898 4.54178 2.38735
C 32.57095 5.49162 3.34597
C 33.45908 6.55380 2.84881
C 33.84018 6.22888 1.44741
C 32.98687 5.08381 1.20092
C 33.92667 7.69216 3.70841
C 34.80276 6.85410 0.48370
O 34.79970 6.55266 -0.76208
C 35.68626 7.89779 1.02373
N 31.23942 2.72217 0.34897
C 31.87083 3.56115 -0.49088
C 31.71787 3.11456 -1.92530
C 31.05327 1.66168 -1.85439
C 30.84608 1.55242 -0.32778
C 30.86799 4.01972 -2.81681
C 31.79877 0.51041 -2.58341
C 33.10989 0.22274 -1.88800
N 30.39035 1.01667 2.56432
C 29.88878 0.22793 1.50582
C 29.07339 -0.84510 2.11095
C 29.13883 -0.59060 3.48915
C 29.92193 0.57158 3.72490
C 28.35873 -1.84905 1.39860
C 28.71475 -1.03469 4.79745
O 28.10674 -2.01401 5.22755
C 29.05326 0.10509 5.85230
C 29.72586 -0.53975 7.03311
O 30.84246 -1.07345 6.91075
O 29.04406 -0.23118 8.15128
C 29.67895 -0.83894 9.35102
H 31.74634 6.45025 5.04279
H 33.06258 5.05535 -0.90139
H 29.76291 -0.26771 -0.44098
H 29.35848 2.82199 7.38974
H 30.99286 5.10095 7.12098
H 28.22353 4.47017 6.18661
H 29.09936 5.64365 5.18154
H 28.86426 6.01509 6.86273
H 32.36465 3.17517 7.45297
H 31.69833 1.50630 7.57975
H 30.34916 3.14189 9.48031
H 31.95478 3.77578 9.40917
H 33.35037 7.75612 4.63964
H 33.68335 8.61705 3.17434
H 35.00227 7.66456 3.89544
H 35.08345 8.70141 1.43998
H 36.36169 8.16020 0.20714
H 36.26193 7.47369 1.85401
H 32.68281 3.03628 -2.41095
H 30.05183 1.64908 -2.30913
H 31.40689 4.37951 -3.71359
H 30.51920 4.85387 -2.22261
H 29.93579 3.57959 -3.16680
H 32.04799 0.80774 -3.60500
H 31.22169 -0.40610 -2.64646
H 33.78315 0.97507 -2.30968
H 33.45151 -0.78487 -2.11615
H 33.03460 0.41422 -0.81128
H 28.53157 -2.78216 1.91815
H 28.68053 -1.90545 0.36079
H 27.33904 -1.46594 1.49193
H 28.11380 0.58552 6.12829
H 28.90985 -1.43471 9.85225
H 30.11237 0.00521 9.91961
H 30.44970 -1.58301 9.18750

