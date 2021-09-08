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
Mg 7.89043 0.70477 0.60473
C 8.09537 -2.46005 -1.07484
C 6.94915 -1.00838 3.40003
C 7.26570 3.49366 1.98287
C 7.72496 2.20438 -2.53684
N 7.59663 -1.45362 0.98851
C 7.76381 -2.56767 0.29188
C 7.64314 -3.87164 1.12534
C 6.73626 -3.30747 2.32457
C 7.13960 -1.83561 2.26445
C 5.24504 -3.40571 2.04898
C 9.03359 -4.37889 1.64628
C 9.03101 -5.80479 2.20176
H 10.21132 -6.63746 1.98729
N 7.68639 1.17245 2.52260
C 7.31900 0.34866 3.52643
C 7.29768 1.17879 4.74075
C 7.52224 2.59556 4.34444
C 7.47847 2.47559 2.90091
C 7.15866 0.60084 6.11918
C 7.72793 3.85265 5.13413
O 7.66629 5.00971 4.58619
C 7.91082 3.69342 6.58419
N 7.33275 2.55858 -0.15336
C 7.17589 3.61876 0.65863
C 7.04601 4.89794 -0.13335
C 7.39194 4.50470 -1.64440
C 7.58331 2.98185 -1.47211
C 5.67615 5.57434 -0.08463
C 8.51681 5.31666 -2.34303
C 9.84864 5.01817 -1.69301
N 8.08957 0.04895 -1.44424
C 7.95840 0.83700 -2.60866
C 8.06378 -0.07522 -3.76588
C 8.22721 -1.34391 -3.18950
C 8.20368 -1.23300 -1.77294
C 7.95749 0.30143 -5.13444
C 8.40364 -2.75311 -3.45729
O 8.58363 -3.42915 -4.46938
C 8.19942 -3.55485 -2.10013
C 9.33171 -4.53578 -1.96523
O 10.50011 -4.13603 -1.81807
O 8.82664 -5.77372 -1.81435
C 9.91107 -6.78053 -1.66550
H 6.42909 -1.42420 4.26347
H 7.15548 4.49476 2.40609
H 7.78137 2.60712 -3.55040
H 7.13275 -4.68542 0.62430
H 7.02543 -3.74346 3.27728
H 5.04248 -3.73118 1.03943
H 4.84955 -2.41179 2.21260
H 4.76551 -4.04727 2.79381
H 9.29973 -3.82374 2.55497
H 9.78201 -4.28978 0.88302
H 8.12020 -6.35168 1.92919
H 8.99571 -5.75217 3.29263
H 6.84492 -0.44978 6.08510
H 6.32182 1.11874 6.60016
H 8.05581 0.73914 6.72636
H 7.04325 3.19032 7.00452
H 8.13222 4.69334 6.96236
H 8.77048 3.03671 6.75836
H 7.75826 5.63867 0.20887
H 6.52589 4.63940 -2.30896
H 5.71772 6.61872 0.27814
H 5.02407 4.97997 0.54150
H 5.13596 5.59199 -1.02966
H 8.34375 6.38838 -2.21811
H 8.59052 5.12434 -3.40819
H 9.85508 5.68537 -0.82590
H 10.66667 5.24691 -2.37312
H 9.88129 3.98891 -1.31722
H 8.74101 -0.22427 -5.66380
H 8.04650 1.37869 -5.25839
H 6.94834 -0.05538 -5.35712
H 7.23077 -4.05296 -2.15647
H 9.75230 -7.53789 -2.43937
H 9.86304 -7.11796 -0.61308
H 10.92360 -6.46275 -1.88515

