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
Mg 8.94058 0.81221 0.68182
C 8.47455 0.68258 4.23781
C 7.76092 4.01072 0.79399
C 8.91302 0.80782 -2.49098
C 8.98380 -2.67254 0.70009
N 8.24069 2.07319 2.35882
C 8.16685 1.93850 3.67441
C 7.81663 3.24495 4.43605
C 7.07183 4.02475 3.24607
C 7.75373 3.35851 2.05283
C 5.58310 3.72891 3.17184
C 9.09490 4.01727 4.91726
C 8.82874 5.13355 5.92960
H 9.81665 5.37233 6.97817
N 8.88601 2.32284 -0.60411
C 8.39275 3.55989 -0.38579
C 8.56852 4.29194 -1.64967
C 9.05495 3.33625 -2.68154
C 8.94801 2.08692 -1.95504
C 8.35778 5.77350 -1.76797
C 9.52491 3.51539 -4.09336
O 9.67569 2.51828 -4.88449
C 9.71374 4.89980 -4.55079
N 8.73564 -0.70539 -0.72433
C 8.81229 -0.44300 -2.04084
C 8.91322 -1.71659 -2.84591
C 9.13288 -2.88485 -1.77610
C 9.02709 -2.05959 -0.47489
C 7.70427 -2.04838 -3.72040
C 10.37490 -3.79631 -1.97438
C 11.64021 -3.00171 -1.74420
N 8.94980 -0.73528 2.18954
C 8.94284 -2.13157 1.97879
C 8.83501 -2.77164 3.30585
C 8.76296 -1.69450 4.20199
C 8.80140 -0.46958 3.48254
C 8.76614 -4.17213 3.55083
C 8.65243 -1.31519 5.59217
O 8.66806 -1.92233 6.66216
C 8.33336 0.24008 5.66780
C 9.25498 0.85436 6.68554
O 10.48270 0.88468 6.48960
O 8.52170 1.49491 7.61418
C 9.39173 2.13003 8.63961
H 7.19397 4.93576 0.68569
H 9.01290 0.74886 -3.57727
H 9.08949 -3.75509 0.79802
H 7.14289 3.11525 5.27452
H 7.29666 5.08790 3.26900
H 5.29246 2.94923 3.86015
H 5.39465 3.42223 2.15131
H 5.00833 4.64736 3.32154
H 9.49001 4.61074 4.08274
H 9.82454 3.33750 5.31243
H 7.82119 5.06968 6.35821
H 8.83508 6.08998 5.40120
H 7.84355 6.17823 -0.88754
H 7.65241 5.93143 -2.59099
H 9.28134 6.31427 -1.98538
H 8.77660 5.44292 -4.45457
H 10.13719 4.82093 -5.55396
H 10.43270 5.39265 -3.88685
H 9.76607 -1.68594 -3.51287
H 8.29166 -3.59345 -1.76808
H 7.95971 -2.17124 -4.78985
H 6.96626 -1.26758 -3.59298
H 7.15166 -2.93890 -3.42576
H 10.41863 -4.15798 -3.00464
H 10.37957 -4.66686 -1.32700
H 11.80153 -2.50986 -2.70812
H 12.46776 -3.66287 -1.49505
H 11.48121 -2.21650 -0.99599
H 9.41637 -4.37345 4.39176
H 9.06086 -4.74417 2.67351
H 7.70081 -4.27930 3.77188
H 7.28473 0.34611 5.94835
H 9.06605 1.75181 9.61356
H 9.30813 3.21980 8.46885
H 10.44031 1.85638 8.64633

