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
Mg 8.94228 0.79783 0.67640
C 8.98832 1.00370 -2.90613
C 8.38693 -2.54951 0.32784
C 8.44837 0.65455 3.80737
C 8.35996 4.23352 0.72774
N 8.73493 -0.48948 -1.11024
C 8.83223 -0.31076 -2.41903
C 8.83861 -1.62273 -3.24863
C 8.07961 -2.57368 -2.20054
C 8.44760 -1.85268 -0.90538
C 6.56748 -2.54681 -2.34911
C 10.29160 -2.13678 -3.54292
C 10.38238 -3.23572 -4.60398
H 11.54033 -3.25034 -5.49358
N 8.97203 -0.75531 1.91149
C 8.74696 -2.04854 1.59817
C 8.86388 -2.79460 2.86067
C 9.01337 -1.81591 3.97176
C 8.79133 -0.57321 3.26008
C 8.90812 -4.29397 2.91870
C 9.29522 -1.97331 5.43531
O 9.14617 -1.00183 6.25802
C 9.66075 -3.32240 5.89117
N 8.26321 2.19124 2.06190
C 8.19148 1.88736 3.36974
C 7.94216 3.12141 4.20361
C 8.10283 4.35702 3.20109
C 8.33967 3.58538 1.88432
C 6.57806 3.19663 4.88914
C 9.11771 5.46207 3.60304
C 10.52521 4.91255 3.55442
N 8.89440 2.38926 -0.78411
C 8.60540 3.75183 -0.55175
C 8.58039 4.42218 -1.86797
C 8.83584 3.39050 -2.78387
C 8.98790 2.16014 -2.08918
C 8.29740 5.79853 -2.09531
C 9.00091 3.06063 -4.18124
O 9.06412 3.70865 -5.22519
C 8.98154 1.47852 -4.33271
C 10.13805 1.08159 -5.20874
O 11.30884 1.25799 -4.82827
O 9.67659 0.36496 -6.24983
C 10.78780 -0.06116 -7.14171
H 7.98603 -3.56349 0.33164
H 8.37526 0.68115 4.89706
H 8.28236 5.32114 0.66662
H 8.28321 -1.57545 -4.17765
H 8.49296 -3.57888 -2.20804
H 6.24566 -1.79992 -3.05955
H 6.17901 -2.32411 -1.36394
H 6.19571 -3.54406 -2.60132
H 10.66003 -2.68880 -2.66869
H 10.93703 -1.32263 -3.80976
H 9.45390 -3.33009 -5.18014
H 10.48298 -4.19905 -4.09812
H 8.61018 -4.74248 1.96294
H 8.12950 -4.60988 3.62142
H 9.87176 -4.67385 3.26476
H 8.86113 -4.01638 5.64326
H 9.91099 -3.21580 6.94848
H 10.54631 -3.65131 5.33599
H 8.67916 3.21062 4.99227
H 7.15830 4.90695 3.07753
H 6.64723 3.31408 5.98709
H 6.01947 2.30517 4.63651
H 5.92364 3.98893 4.52970
H 8.94374 5.77907 4.63425
H 9.06075 6.34822 2.97975
H 10.62894 4.41369 4.52265
H 11.24776 5.71897 3.44647
H 10.62182 4.14611 2.77668
H 9.01721 6.14826 -2.82337
H 8.35219 6.37329 -1.17313
H 7.27443 5.72729 -2.47443
H 8.02192 1.20319 -4.77191
H 10.54613 0.29760 -8.14698
H 10.87749 -1.15529 -7.00538
H 11.75944 0.39194 -6.98336

