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
Mg 15.77583 1.42391 1.19860
C 17.28556 -0.69494 3.67050
C 16.23046 4.02076 3.36292
C 14.05107 3.22838 -0.76010
C 14.43159 -1.45016 -0.24308
N 16.65645 1.52427 3.22510
C 17.25386 0.66635 4.03837
C 17.92161 1.31556 5.28022
C 17.03899 2.65506 5.35513
C 16.64115 2.76817 3.88489
C 15.77357 2.49654 6.18152
C 19.43833 1.63962 5.04240
C 20.22752 1.98744 6.30650
H 21.62036 1.55829 6.39761
N 15.64611 3.39880 1.05183
C 15.88827 4.30456 2.02251
C 15.63209 5.62083 1.41755
C 15.06023 5.40980 0.05992
C 14.87350 3.97277 0.07298
C 16.00113 6.91117 2.09004
C 14.72366 6.34783 -1.05947
O 14.01472 5.97516 -2.06014
C 15.16944 7.74022 -0.90445
N 14.26510 0.98411 -0.16030
C 13.72009 1.94738 -0.92384
C 12.85658 1.35696 -2.01294
C 13.12938 -0.21829 -1.97096
C 14.07715 -0.27861 -0.75297
C 11.35565 1.61682 -1.88616
C 13.58472 -0.88924 -3.29579
C 14.96596 -0.40353 -3.67215
N 15.97615 -0.70633 1.49865
C 15.30245 -1.73301 0.80143
C 15.67968 -3.01165 1.43798
C 16.53540 -2.63528 2.48420
C 16.66521 -1.22039 2.51127
C 15.20277 -4.29962 1.06408
C 17.34999 -3.11245 3.57849
O 17.69425 -4.22110 3.98606
C 17.77179 -1.86896 4.47414
C 19.24813 -1.97024 4.74377
O 20.06872 -1.86852 3.81484
O 19.45553 -1.95228 6.07306
C 20.90589 -2.03935 6.39042
H 16.11840 4.85837 4.05193
H 13.51905 3.79590 -1.52715
H 14.13217 -2.39817 -0.69509
H 17.81854 0.75052 6.19885
H 17.64152 3.50891 5.65387
H 15.61647 1.47226 6.48535
H 14.96372 2.83129 5.54658
H 15.79151 3.18199 7.03367
H 19.52064 2.59126 4.50171
H 19.92296 0.83546 4.52366
H 19.68228 1.71939 7.21963
H 20.32916 3.07388 6.36419
H 16.25187 6.75668 3.14676
H 15.09739 7.52971 2.11294
H 16.78690 7.45304 1.55936
H 14.73970 8.15811 0.00283
H 14.91490 8.23486 -1.84382
H 16.25680 7.74647 -0.76949
H 13.14557 1.73787 -2.98479
H 12.22577 -0.78106 -1.69417
H 10.92302 2.11784 -2.77266
H 11.18539 2.20652 -0.99521
H 10.74826 0.73317 -1.69790
H 12.92056 -0.60070 -4.11420
H 13.58928 -1.97322 -3.25166
H 14.76443 0.54498 -4.17890
H 15.45251 -1.11199 -4.33952
H 15.56147 -0.18023 -2.77931
H 16.05783 -4.96244 1.07141
H 14.72379 -4.27813 0.08741
H 14.48595 -4.48338 1.86886
H 17.18412 -1.91195 5.39195
H 21.03568 -2.89557 7.05954
H 21.18555 -1.04455 6.78529
H 21.58099 -2.29768 5.58299

