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
Mg 13.14655 1.17592 0.99719
C 13.25161 0.70374 -2.55880
C 12.54361 -2.16993 1.27347
C 12.60175 1.63033 4.08978
C 12.61814 4.56706 0.39164
N 12.94662 -0.42171 -0.51935
C 13.06711 -0.49357 -1.83641
C 13.06557 -1.93813 -2.40456
C 12.27533 -2.66529 -1.21056
C 12.63452 -1.71829 -0.06737
C 10.76632 -2.64768 -1.38913
C 14.51461 -2.51672 -2.57045
C 14.60443 -3.79666 -3.40437
H 15.77571 -3.99302 -4.25409
N 13.13241 -0.11746 2.50235
C 12.89172 -1.44359 2.43346
C 12.97711 -1.94031 3.81554
C 13.12482 -0.77209 4.72557
C 12.93365 0.31730 3.78920
C 12.99662 -3.40245 4.15496
C 13.38131 -0.65498 6.19754
O 13.23492 0.45570 6.82027
C 13.71823 -1.89879 6.90529
N 12.46821 2.81353 2.08366
C 12.37132 2.76198 3.42376
C 12.12868 4.13390 4.00630
C 12.32453 5.15677 2.79265
C 12.56956 4.14832 1.64881
C 10.75543 4.35398 4.64045
C 13.35045 6.30475 2.99836
C 14.74964 5.73810 3.07955
N 13.14669 2.46475 -0.73687
C 12.87579 3.85025 -0.76991
C 12.88192 4.26135 -2.18884
C 13.13518 3.07267 -2.88984
C 13.25685 1.99310 -1.97375
C 12.62440 5.57387 -2.67579
C 13.31672 2.48382 -4.19712
O 13.40646 2.92308 -5.34287
C 13.27459 0.90183 -4.04906
C 14.43815 0.33260 -4.81366
O 15.60552 0.56256 -4.45176
O 13.98162 -0.56114 -5.70986
C 15.09968 -1.16144 -6.48537
H 12.12664 -3.15993 1.46033
H 12.51211 1.86232 5.15356
H 12.55879 5.62467 0.12591
H 12.52551 -2.05938 -3.33591
H 12.67274 -3.65912 -1.02156
H 10.46750 -2.04372 -2.23298
H 10.36615 -2.23877 -0.47060
H 10.38272 -3.66977 -1.45630
H 14.86058 -2.89910 -1.60152
H 15.17697 -1.77553 -2.97363
H 13.68365 -3.98594 -3.96943
H 14.68184 -4.64886 -2.72483
H 12.70652 -4.01889 3.29517
H 12.20224 -3.57066 4.89016
H 13.94860 -3.72265 4.58377
H 12.91164 -2.61682 6.77757
H 13.95364 -1.59842 7.92815
H 14.60699 -2.33744 6.43805
H 12.85463 4.36050 4.77750
H 11.39089 5.68558 2.55075
H 10.80936 4.67495 5.69788
H 10.18676 3.43807 4.54960
H 10.11935 5.07278 4.12670
H 13.16551 6.81222 3.94832
H 13.31729 7.05851 2.21879
H 14.83035 5.42897 4.12598
H 15.48650 6.50060 2.83528
H 14.84616 4.83790 2.46151
H 13.36092 5.77128 -3.44330
H 12.67395 6.31106 -1.87713
H 11.60646 5.44556 -3.05341
H 12.31767 0.56100 -4.44618
H 14.87941 -0.99511 -7.54438
H 15.16984 -2.21145 -6.14434
H 16.07580 -0.69899 -6.39722

