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
Mg 8.92994 0.80857 0.68739
C 6.64683 -1.78689 1.65178
C 6.23207 2.87945 0.42752
C 10.84948 2.98543 -0.59482
C 11.23309 -1.69279 -0.07709
N 6.76957 0.47335 1.02308
C 6.02574 -0.53715 1.44727
C 4.54696 -0.16814 1.74156
C 4.42653 1.11977 0.78990
C 5.89313 1.54558 0.76812
C 4.00634 0.78117 -0.63071
C 4.31259 0.20626 3.24719
C 2.84393 0.27671 3.67125
H 2.48729 -0.15586 5.01964
N 8.69260 2.76804 0.48077
C 7.52640 3.44418 0.41666
C 7.89057 4.85880 0.24237
C 9.36302 4.94306 0.04358
C 9.69675 3.53623 -0.05410
C 6.89328 5.97533 0.35274
C 10.31174 6.09896 -0.05651
O 11.51626 5.94777 -0.46774
C 9.76177 7.42708 0.25188
N 10.69342 0.66557 -0.40442
C 11.35495 1.76845 -0.79678
C 12.72795 1.42303 -1.32208
C 12.94440 -0.12297 -0.97481
C 11.56291 -0.44040 -0.36154
C 12.94200 1.64282 -2.81950
C 14.21601 -0.49036 -0.16196
C 14.11384 0.07116 1.23791
N 9.03662 -1.33534 0.93303
C 10.08126 -2.19043 0.51877
C 9.64823 -3.57460 0.79981
C 8.36207 -3.42898 1.34117
C 8.00975 -2.05268 1.37512
C 10.39113 -4.75370 0.51009
C 7.21809 -4.12736 1.88180
O 6.97440 -5.29629 2.17835
C 6.00630 -3.10448 1.98990
C 5.39169 -3.25502 3.35449
O 6.03527 -2.95173 4.37450
O 4.08010 -3.52328 3.21986
C 3.42421 -3.67037 4.54651
H 5.43725 3.54823 0.09605
H 11.59300 3.70268 -0.95007
H 11.93340 -2.52229 -0.19608
H 3.82387 -0.91847 1.44504
H 3.81209 1.89174 1.24586
H 3.97251 -0.28521 -0.79750
H 4.74518 1.24142 -1.27362
H 3.05555 1.26676 -0.86846
H 4.60428 1.25248 3.40559
H 4.84362 -0.46495 3.89381
H 2.17864 -0.18438 2.93117
H 2.53807 1.32558 3.68993
H 5.86344 5.59767 0.33877
H 6.97681 6.57039 -0.56311
H 7.07884 6.62236 1.21274
H 8.93001 7.63754 -0.41621
H 10.61136 8.11156 0.21234
H 9.35237 7.40886 1.26808
H 13.49142 2.00587 -0.82155
H 13.03318 -0.73617 -1.88362
H 13.78635 2.32136 -3.04475
H 12.02304 2.02572 -3.24308
H 13.10434 0.73670 -3.40070
H 15.09993 -0.03502 -0.61527
H 14.39786 -1.55832 -0.10377
H 14.44744 1.10575 1.11375
H 14.76478 -0.47453 1.91800
H 13.07132 0.09925 1.57557
H 10.30083 -5.39872 1.37397
H 11.43193 -4.52416 0.29141
H 9.85492 -5.10986 -0.37359
H 5.30091 -3.34700 1.19417
H 2.93679 -4.65018 4.55506
H 2.77322 -2.78293 4.65796
H 4.06534 -3.74549 5.41700

