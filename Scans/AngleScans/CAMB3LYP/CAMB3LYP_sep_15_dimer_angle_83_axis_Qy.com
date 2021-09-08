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
Mg 7.89095 0.71327 0.59262
C 10.56859 -1.47887 1.54333
C 9.84503 3.30798 1.63362
C 5.39254 2.62692 0.18849
C 5.77235 -2.05154 0.70681
N 9.92561 0.77907 1.45754
C 10.87101 -0.11200 1.71555
C 12.24233 0.50047 2.10786
C 11.72521 1.92142 2.64856
C 10.43392 2.03505 1.84089
C 11.36557 1.90887 4.12499
C 13.19699 0.66825 0.87397
C 14.65552 0.96955 1.22561
H 15.70557 0.40859 0.37992
N 7.80209 2.68694 0.40474
C 8.69235 3.58933 0.86778
C 8.16203 4.90834 0.48945
C 6.80435 4.71945 -0.08998
C 6.60075 3.30841 0.17007
C 8.96040 6.17323 0.61625
C 5.84249 5.65680 -0.75511
O 4.62159 5.32969 -0.96775
C 6.34698 7.00238 -1.06558
N 5.83894 0.38866 0.67145
C 4.97134 1.38069 0.40487
C 3.56594 0.84980 0.25191
C 3.70816 -0.74298 0.23182
C 5.22821 -0.86442 0.47695
C 2.57787 1.25853 1.34418
C 3.09358 -1.48915 -0.98398
C 3.86357 -1.15085 -2.24019
N 8.12787 -1.41791 0.85719
C 7.10564 -2.39131 0.89730
C 7.74747 -3.68393 1.21292
C 9.10632 -3.36634 1.35838
C 9.29498 -1.97073 1.16817
C 7.07583 -4.92822 1.37649
C 10.42251 -3.89264 1.64007
O 10.89244 -5.01881 1.79707
C 11.40975 -2.67379 1.89721
C 12.66044 -2.91214 1.09628
O 12.62395 -2.91519 -0.14681
O 13.72390 -2.88036 1.91994
C 14.98925 -3.10080 1.17021
H 10.28138 4.17000 2.13900
H 4.51073 3.22387 -0.05552
H 5.19413 -2.97699 0.66350
H 12.76750 -0.02795 2.89447
H 12.41314 2.71958 2.38203
H 11.40509 0.91305 4.54083
H 10.36033 2.30484 4.18627
H 12.00001 2.61097 4.67350
H 12.93681 1.59303 0.34309
H 13.14957 -0.19258 0.23572
H 14.87204 0.77976 2.28394
H 14.82713 2.04195 1.10506
H 9.85904 6.02219 1.22701
H 8.35359 6.87766 1.19541
H 9.19548 6.62173 -0.35136
H 6.68007 7.48354 -0.14907
H 5.54435 7.49522 -1.61764
H 7.22884 6.90562 -1.70868
H 3.12928 1.17546 -0.68424
H 3.21226 -1.20727 1.09690
H 1.68265 1.77653 0.95133
H 3.09725 1.88481 2.05732
H 2.21895 0.44373 1.97066
H 2.06569 -1.15814 -1.15138
H 3.06901 -2.56688 -0.86195
H 3.42056 -0.20052 -2.55278
H 3.72070 -1.92048 -2.99603
H 4.92018 -0.96402 -2.01586
H 7.66568 -5.67115 0.85627
H 6.05938 -4.88508 0.99081
H 7.09803 -5.02014 2.46568
H 11.61076 -2.63540 2.96851
H 15.49656 -3.94811 1.64179
H 15.51644 -2.12833 1.18141
H 14.91124 -3.44468 0.14546

