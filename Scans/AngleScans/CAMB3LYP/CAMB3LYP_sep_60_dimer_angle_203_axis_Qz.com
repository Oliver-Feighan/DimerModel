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
Mg 31.54045 2.83812 2.38829
C 31.18256 5.96417 0.66237
C 31.02393 1.29478 -0.60936
C 31.40168 0.07743 3.94610
C 30.89160 4.55755 5.34941
N 31.16390 3.63458 0.35947
C 31.11082 4.84011 -0.18663
C 31.05378 4.84157 -1.73785
C 30.43748 3.37497 -1.95720
C 30.93141 2.70802 -0.67506
C 28.91799 3.35372 -1.95185
C 32.47595 4.97016 -2.38816
C 32.47226 5.28828 -3.88500
H 33.52395 6.15125 -4.41576
N 31.74196 0.97821 1.72594
C 31.51205 0.53102 0.47357
C 31.78536 -0.91411 0.50587
C 32.03600 -1.31266 1.91763
C 31.70862 -0.07803 2.60217
C 31.86895 -1.75438 -0.73525
C 32.48254 -2.59675 2.54843
O 32.40461 -2.78836 3.81336
C 32.92634 -3.65759 1.63238
N 30.98928 2.34119 4.33035
C 31.07043 1.07272 4.76890
C 30.88831 1.00052 6.26623
C 30.90892 2.51650 6.77512
C 31.01463 3.23889 5.41411
C 29.60773 0.31943 6.74853
C 31.93756 2.87303 7.88303
C 33.34270 2.77143 7.33471
N 31.29497 4.90595 2.96436
C 31.00358 5.40348 4.25340
C 30.81728 6.86333 4.12477
C 30.99188 7.11260 2.75515
C 31.24756 5.89138 2.07490
C 30.48035 7.75529 5.18166
C 31.00751 8.12973 1.72854
O 30.93969 9.35818 1.72170
C 31.00640 7.41490 0.30889
C 32.06663 8.06330 -0.53847
O 33.26832 7.96131 -0.23490
O 31.50803 8.50906 -1.67849
C 32.52213 9.15047 -2.55713
H 30.64847 0.71064 -1.45001
H 31.44984 -0.83231 4.54904
H 30.78364 5.17823 6.24147
H 30.39743 5.58878 -2.16758
H 30.86956 2.89135 -2.82943
H 28.50259 4.31345 -1.68260
H 28.63685 2.60253 -1.22540
H 28.54259 2.99677 -2.91518
H 32.95160 3.98116 -2.40744
H 33.06947 5.69528 -1.86609
H 31.48769 5.62406 -4.23261
H 32.64988 4.36215 -4.43704
H 31.47640 -1.22064 -1.60957
H 31.18063 -2.59551 -0.59890
H 32.87331 -2.14480 -0.91274
H 32.11962 -3.90251 0.94558
H 33.29037 -4.46332 2.27282
H 33.75160 -3.27156 1.02366
H 31.70642 0.46305 6.72984
H 29.94457 2.81219 7.21378
H 29.79616 -0.53392 7.42704
H 29.04434 -0.00081 5.88218
H 28.89998 0.97290 7.25574
H 31.87258 2.15597 8.70507
H 31.79194 3.86154 8.30569
H 33.56466 1.70354 7.42101
H 34.03071 3.36704 7.93130
H 33.36900 3.02854 6.26941
H 31.10666 8.62943 5.06291
H 30.62479 7.29102 6.15499
H 29.42334 7.92766 4.96190
H 30.01016 7.53613 -0.11816
H 32.16210 10.16047 -2.77614
H 32.65063 8.46438 -3.41533
H 33.49510 9.36102 -2.12880

