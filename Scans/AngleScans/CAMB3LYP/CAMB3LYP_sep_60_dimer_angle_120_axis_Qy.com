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
Mg 31.53538 2.84157 2.38549
C 34.14654 0.46671 1.73668
C 33.85314 5.29703 1.90249
C 29.37697 4.93909 3.39002
C 29.75702 0.26079 3.90954
N 33.68920 2.76606 1.88788
C 34.55957 1.80882 1.60418
C 35.92704 2.31890 1.07558
C 35.89883 3.79091 1.71690
C 34.38718 3.98673 1.81299
C 36.47791 3.84692 3.12078
C 35.97844 2.37603 -0.49161
C 37.37688 2.56795 -1.08262
H 37.69781 1.90165 -2.34179
N 31.45074 4.81023 2.14938
C 32.48623 5.64990 1.93960
C 31.90093 6.99681 1.85346
C 30.45521 6.90266 2.19367
C 30.37433 5.52013 2.62046
C 32.68217 8.19608 1.40088
C 29.33401 7.89645 2.15502
O 28.20783 7.66336 2.72091
C 29.62462 9.18814 1.51583
N 29.91071 2.68583 3.67333
C 29.10302 3.73716 3.89757
C 27.85311 3.31752 4.63379
C 27.87713 1.71849 4.64437
C 29.24151 1.48137 3.96053
C 27.71989 3.83372 6.06633
C 26.62828 0.99209 4.07396
C 26.52602 1.23342 2.58510
N 31.77681 0.70606 2.60771
C 30.92777 -0.18025 3.30605
C 31.56737 -1.51154 3.27440
C 32.76476 -1.30087 2.57413
C 32.87396 0.06937 2.21357
C 31.06042 -2.69275 3.88585
C 33.96621 -1.92395 2.06719
O 34.38196 -3.07977 1.99661
C 34.97428 -0.78230 1.61244
C 35.49939 -1.14214 0.24962
O 34.73751 -1.17549 -0.73274
O 36.84391 -1.17427 0.28982
C 37.41189 -1.51609 -1.04149
H 34.54547 6.13423 1.99553
H 28.55156 5.59851 3.66811
H 29.21933 -0.61581 4.27716
H 36.78784 1.77107 1.43983
H 36.33619 4.52354 1.04360
H 36.70545 2.86264 3.50230
H 35.72284 4.32422 3.73161
H 37.34750 4.50996 3.14368
H 35.50164 3.30426 -0.83191
H 35.52158 1.50507 -0.91980
H 38.16557 2.38968 -0.34158
H 37.49736 3.61939 -1.35460
H 33.75928 7.98941 1.37831
H 32.56881 8.96242 2.17531
H 32.32396 8.59810 0.45076
H 30.45701 9.66520 2.02769
H 28.67644 9.72920 1.50597
H 29.95223 9.00336 0.48669
H 26.96557 3.65219 4.11096
H 27.96380 1.31904 5.66553
H 26.79207 4.41158 6.23754
H 28.58990 4.43490 6.29482
H 27.75911 3.06760 6.83881
H 25.71701 1.40037 4.51786
H 26.62700 -0.07649 4.26139
H 26.03163 2.20769 2.52641
H 25.92742 0.45788 2.11154
H 27.51963 1.34031 2.13457
H 31.19291 -3.49451 3.17162
H 30.01555 -2.57806 4.16666
H 31.71542 -2.75699 4.75884
H 35.76941 -0.73161 2.35727
H 38.05691 -2.38865 -0.89917
H 37.89183 -0.58938 -1.40840
H 36.72826 -1.87984 -1.79964

