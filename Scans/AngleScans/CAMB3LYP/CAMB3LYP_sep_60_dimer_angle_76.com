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
Mg 31.53979 2.83422 2.40157
C 31.73829 -0.43104 0.92586
C 30.57460 1.30493 5.29368
C 30.92131 5.70724 3.59768
C 31.39911 4.13367 -0.82912
N 31.23316 0.70567 2.91879
C 31.39866 -0.45081 2.29469
C 31.26692 -1.69912 3.20789
C 30.35624 -1.05634 4.36391
C 30.76721 0.40694 4.21358
C 28.86610 -1.16487 4.08638
C 32.65194 -2.17893 3.76783
C 32.63921 -3.56698 4.41197
H 33.81654 -4.41693 4.25725
N 31.32749 3.42274 4.28496
C 30.95048 2.66553 5.33651
C 30.92658 3.57060 6.49601
C 31.16035 4.95853 6.01258
C 31.12396 4.74806 4.57924
C 30.77709 3.08129 7.90727
C 31.36792 6.26192 6.72272
O 31.31505 7.38244 6.10265
C 31.54201 6.19352 8.18097
N 30.99551 4.63914 1.52500
C 30.83943 5.74908 2.26769
C 30.72028 6.97640 1.39599
C 31.07260 6.48716 -0.08526
C 31.25544 4.97732 0.18370
C 29.35355 7.66083 1.39403
C 32.20532 7.24829 -0.82708
C 33.53204 6.98522 -0.15181
N 31.74698 2.04971 0.39916
C 31.62615 2.76342 -0.81332
C 31.73339 1.77963 -1.91017
C 31.88734 0.54904 -1.25410
C 31.85654 0.74908 0.15251
C 31.63652 2.06980 -3.30033
C 32.05824 -0.87504 -1.43158
O 32.24046 -1.61431 -2.39802
C 31.84256 -1.58873 -0.02783
C 32.96921 -2.56441 0.17517
O 34.13875 -2.16157 0.30366
O 32.45717 -3.78805 0.40076
C 33.53575 -4.78847 0.61903
H 30.04772 0.94674 6.17854
H 30.81373 6.73351 3.95637
H 31.46311 4.47149 -1.86568
H 30.75527 -2.54048 2.75613
H 30.63798 -1.43277 5.34385
H 28.66751 -1.55236 3.09817
H 28.47465 -0.16081 4.18478
H 28.37928 -1.75602 4.86733
H 32.91581 -1.56887 4.64129
H 33.40499 -2.14153 3.00484
H 31.72722 -4.12576 4.16910
H 32.59815 -3.44558 5.49714
H 30.45833 2.03206 7.93758
H 29.94018 3.63231 8.34979
H 31.67155 3.25343 8.50974
H 30.66965 5.72190 8.62708
H 31.76628 7.21425 8.49670
H 32.39743 5.54514 8.40115
H 31.43430 7.73393 1.69502
H 30.21090 6.58371 -0.76202
H 29.39830 8.72579 1.69055
H 28.69508 7.11010 2.05255
H 28.81867 7.62139 0.44664
H 32.03689 8.32655 -0.77091
H 32.28394 6.98892 -1.87757
H 33.53701 7.70569 0.67159
H 34.35494 7.16689 -0.84020
H 33.55750 5.98155 0.28824
H 32.42034 1.50819 -3.79096
H 31.73157 3.13669 -3.49135
H 30.62687 1.70432 -3.50596
H 30.87177 -2.08493 -0.05833
H 33.37750 -5.59234 -0.10652
H 33.48024 -5.05870 1.69031
H 34.55105 -4.48983 0.38571

