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
Mg 13.14487 1.19009 1.00252
C 12.74559 0.38746 4.47749
C 11.92703 4.29145 1.73260
C 13.05939 1.78128 -2.11367
C 13.23258 -2.23494 0.36441
N 12.45983 2.10235 2.89745
C 12.41176 1.72180 4.16521
C 12.05898 2.85621 5.16433
C 11.28274 3.83371 4.15398
C 11.95110 3.41437 2.84627
C 9.79677 3.53347 4.04862
C 13.33594 3.54460 5.76224
C 13.07421 4.44649 6.97048
H 14.07801 4.49973 8.02974
N 13.04768 2.91430 0.02466
C 12.54285 4.08034 0.47939
C 12.68620 5.03944 -0.62683
C 13.16576 4.30244 -1.82748
C 13.08796 2.93734 -1.34733
C 12.45458 6.51332 -0.46109
C 13.60754 4.75101 -3.18756
O 13.75646 3.92277 -4.15434
C 13.77042 6.19952 -3.37933
N 12.93352 -0.03945 -0.66059
C 12.98277 0.46674 -1.90528
C 13.08511 -0.63120 -2.93700
C 13.33907 -1.97597 -2.10954
C 13.24662 -1.41159 -0.67486
C 11.86467 -0.81199 -3.83935
C 14.58870 -2.81412 -2.49502
C 15.84784 -2.05691 -2.13923
N 13.20126 -0.61271 2.19197
C 13.20814 -1.94449 1.72250
C 13.13270 -2.82407 2.90702
C 13.06340 -1.93570 3.99076
C 13.07316 -0.59700 3.51401
C 13.08607 -4.24656 2.88527
C 12.97349 -1.82606 5.42902
O 13.01636 -2.62304 6.36535
C 12.63617 -0.31794 5.80076
C 13.56838 0.10880 6.90140
O 14.79183 0.19488 6.69557
O 12.84414 0.55179 7.94524
C 13.72464 0.99673 9.05815
H 11.34651 5.21125 1.80905
H 13.14015 1.92901 -3.19308
H 13.35376 -3.31480 0.25533
H 11.40238 2.56064 5.97383
H 11.49446 4.87706 4.37295
H 9.52866 2.63387 4.58245
H 9.59360 3.42099 2.99168
H 9.21324 4.39822 4.37732
H 13.70818 4.29045 5.04818
H 14.08124 2.81437 6.01111
H 12.07554 4.28726 7.39506
H 13.05877 5.48513 6.63135
H 11.95144 6.73726 0.48765
H 11.73234 6.81181 -1.22864
H 13.36708 7.09991 -0.58727
H 12.82838 6.69995 -3.16811
H 14.17643 6.31721 -4.38590
H 14.49509 6.57024 -2.64582
H 13.92516 -0.46228 -3.59948
H 12.50717 -2.68673 -2.22184
H 12.10207 -0.72771 -4.91667
H 11.11928 -0.08084 -3.55589
H 11.32887 -1.75066 -3.70889
H 14.61818 -2.97510 -3.57547
H 14.61623 -3.79059 -2.02307
H 15.98527 -1.39026 -2.99584
H 16.68812 -2.73987 -2.03179
H 15.69260 -1.42884 -1.25431
H 13.75406 -4.59189 3.66311
H 13.37194 -4.63888 1.91155
H 12.02641 -4.41025 3.09878
H 11.59159 -0.28317 6.11254
H 13.42164 0.43719 9.94853
H 13.62413 2.09771 9.09671
H 14.77655 0.74338 8.99696

