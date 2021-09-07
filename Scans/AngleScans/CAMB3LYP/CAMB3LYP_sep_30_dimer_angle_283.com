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
Mg 15.77303 1.42729 1.19553
C 15.15748 3.61633 3.97191
C 14.81423 3.98795 -0.84366
C 15.87199 -0.60576 -1.23847
C 15.57400 -1.23821 3.43187
N 15.09575 3.48822 1.62695
C 14.96027 4.22875 2.71672
C 14.67124 5.72911 2.44298
C 14.03037 5.58704 0.97739
C 14.71160 4.29291 0.53713
C 12.52888 5.35362 0.99716
C 15.97969 6.59471 2.41631
C 15.75127 8.10714 2.46236
H 16.71068 8.93415 3.18904
N 15.87464 1.76590 -0.75731
C 15.45993 2.86924 -1.41453
C 15.73632 2.61770 -2.83737
C 16.19621 1.20998 -2.98415
C 15.97406 0.71849 -1.63910
C 15.63361 3.68554 -3.88741
C 16.73347 0.43102 -4.14638
O 16.84630 -0.84514 -4.10720
C 17.03586 1.19616 -5.36482
N 15.51979 -0.63228 1.06752
C 15.66696 -1.27542 -0.10408
C 15.71151 -2.77144 0.09671
C 15.80681 -2.99002 1.67808
C 15.70655 -1.52064 2.14314
C 14.51837 -3.55234 -0.45373
C 16.98963 -3.85126 2.19949
C 18.29674 -3.12885 1.96503
N 15.61482 1.20356 3.33865
C 15.51957 -0.00335 4.06555
C 15.31475 0.35738 5.48342
C 15.28175 1.76001 5.47977
C 15.43365 2.23917 4.15054
C 15.13930 -0.55940 6.55806
C 15.14231 2.94401 6.29672
O 15.07310 3.16208 7.50560
C 14.92898 4.19557 5.34044
C 15.84954 5.29304 5.79938
O 17.08326 5.15690 5.72179
O 15.12589 6.39953 6.04911
C 15.99617 7.51936 6.49681
H 14.31748 4.64472 -1.55836
H 16.01085 -1.34897 -2.02703
H 15.60047 -2.00977 4.20431
H 13.95719 6.18477 3.11847
H 14.32720 6.41189 0.33484
H 12.15767 5.20344 2.00005
H 12.36055 4.47025 0.39519
H 12.01359 6.17066 0.48406
H 16.44798 6.50543 1.42764
H 16.64412 6.30536 3.20716
H 14.72540 8.36031 2.75637
H 15.84495 8.50327 1.44829
H 15.11389 4.57403 -3.50809
H 14.97430 3.29968 -4.67246
H 16.60035 3.93608 -4.32935
H 16.13549 1.70078 -5.70713
H 17.49256 0.48188 -6.05264
H 17.76014 1.97958 -5.11525
H 16.59038 -3.19835 -0.37065
H 14.91887 -3.50571 2.07237
H 14.80720 -4.33816 -1.17701
H 13.83172 -2.85076 -0.90852
H 13.89408 -4.03231 0.29799
H 17.04932 -4.78952 1.64248
H 16.90813 -4.10561 3.25099
H 18.53012 -3.37248 0.92423
H 19.06590 -3.50005 2.63930
H 18.16276 -2.04265 2.02724
H 15.73992 -0.19386 7.38034
H 15.42850 -1.56828 6.27122
H 14.06113 -0.47070 6.71605
H 13.87984 4.48559 5.40956
H 15.60643 7.86123 7.46061
H 15.99514 8.24925 5.66553
H 17.02217 7.28446 6.75487

