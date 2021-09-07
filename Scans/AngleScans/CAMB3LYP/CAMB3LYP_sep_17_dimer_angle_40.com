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
Mg 8.94092 0.80412 0.68689
C 8.99306 -2.70672 1.42874
C 7.77233 1.30351 3.85227
C 8.45632 3.84500 -0.07837
C 9.00533 -0.05018 -2.69123
N 8.48322 -0.59960 2.33381
C 8.60053 -1.90554 2.52120
C 8.35663 -2.37016 3.98209
C 7.44383 -1.14434 4.47506
C 7.95023 -0.06372 3.52191
C 5.96159 -1.35480 4.21430
C 9.68668 -2.46485 4.80914
C 9.56315 -3.20427 6.14321
H 10.69057 -4.01365 6.59740
N 8.69332 2.39625 1.84550
C 8.22994 2.41624 3.11284
C 8.21769 3.83141 3.51446
C 8.55550 4.65952 2.32499
C 8.56122 3.64481 1.29038
C 7.98403 4.27357 4.92985
C 8.81606 6.12452 2.14590
O 8.85666 6.66410 0.98399
C 8.92942 6.92514 3.37376
N 8.54370 1.75800 -1.11718
C 8.42835 3.09611 -1.18096
C 8.41913 3.57518 -2.61302
C 8.79707 2.29662 -3.49621
C 8.87546 1.23215 -2.37990
C 7.09857 4.16374 -3.10910
C 10.00249 2.44186 -4.46500
C 11.28341 2.59171 -3.67623
N 9.17583 -1.01626 -0.45289
C 9.14617 -1.15264 -1.85809
C 9.23425 -2.59664 -2.15732
C 9.28626 -3.20656 -0.89487
C 9.21391 -2.21421 0.11977
C 9.20905 -3.18027 -3.45534
C 9.37521 -4.46497 -0.18975
O 9.54823 -5.63678 -0.52241
C 9.06190 -4.20645 1.34693
C 10.11724 -4.90478 2.16002
O 11.30370 -4.53602 2.10602
O 9.52200 -5.74486 3.02624
C 10.52727 -6.45281 3.86293
H 7.19064 1.55117 4.74063
H 8.39889 4.88766 -0.39936
H 9.12994 -0.39110 -3.72124
H 7.81151 -3.30156 4.07753
H 7.66381 -0.87747 5.50548
H 5.77746 -2.24518 3.63166
H 5.62967 -0.47589 3.67728
H 5.40940 -1.35752 5.15845
H 9.95413 -1.46440 5.17289
H 10.46936 -2.90561 4.22282
H 8.62827 -3.77342 6.21401
H 9.48804 -2.46478 6.94417
H 7.59988 3.45367 5.54919
H 7.16653 5.00228 4.90637
H 8.86364 4.74347 5.37500
H 8.01299 6.83184 3.95161
H 9.20427 7.92887 3.04381
H 9.73389 6.50848 3.99007
H 9.16674 4.34326 -2.76868
H 7.96967 1.99877 -4.15690
H 7.19785 5.19675 -3.49266
H 6.38239 4.12609 -2.29904
H 6.59907 3.58750 -3.88605
H 9.89922 3.34959 -5.06460
H 10.10518 1.61048 -5.15439
H 11.30141 3.65880 -3.43522
H 12.14183 2.30945 -4.28244
H 11.22965 2.04092 -2.72998
H 9.97480 -3.94445 -3.46803
H 9.37735 -2.43481 -4.22977
H 8.18732 -3.56880 -3.47469
H 8.06410 -4.59817 1.54824
H 10.34746 -7.52524 3.73945
H 10.41399 -6.03710 4.88180
H 11.56735 -6.37741 3.56795

