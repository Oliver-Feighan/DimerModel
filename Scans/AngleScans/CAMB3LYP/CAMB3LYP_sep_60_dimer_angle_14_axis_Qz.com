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
Mg 31.53898 2.84196 2.40314
C 31.42384 -0.01489 4.57206
C 30.29382 4.67428 4.99676
C 31.20367 5.27440 0.39361
C 31.67348 0.62426 -0.28189
N 30.97296 2.28837 4.46843
C 31.03067 1.18551 5.19967
C 30.72252 1.39596 6.70643
C 29.84588 2.73622 6.58767
C 30.42566 3.29324 5.28911
C 28.36510 2.46983 6.37403
C 32.02115 1.63229 7.55479
C 31.82587 1.53733 9.06958
H 32.90443 0.97322 9.87640
N 31.32096 4.78198 2.76005
C 30.81954 5.35286 3.87538
C 30.85287 6.80374 3.63454
C 31.26128 7.03556 2.22239
C 31.25758 5.67614 1.72036
C 30.59357 7.81396 4.71419
C 31.58721 8.27672 1.44812
O 31.68609 8.26666 0.17025
C 31.69502 9.52246 2.22143
N 31.23765 2.94221 0.34852
C 31.17936 4.12694 -0.28467
C 31.23442 3.94790 -1.78310
C 31.58690 2.40590 -2.01886
C 31.58671 1.91929 -0.55305
C 29.95523 4.29914 -2.54249
C 32.82686 2.09430 -2.90096
C 34.08766 2.53643 -2.19364
N 31.73437 0.70404 2.15994
C 31.74293 -0.01923 0.94713
C 31.78095 -1.45419 1.29623
C 31.76853 -1.46692 2.69915
C 31.70539 -0.13464 3.18962
C 31.77227 -2.53592 0.37108
C 31.78374 -2.30483 3.87659
O 31.91879 -3.51013 4.08336
C 31.43358 -1.40681 5.14045
C 32.43352 -1.71542 6.22091
O 33.63523 -1.43350 6.06918
O 31.77761 -2.09013 7.33435
C 32.72642 -2.39615 8.43795
H 29.69540 5.29171 5.66724
H 31.19912 6.08081 -0.34339
H 31.81603 -0.12713 -1.06149
H 30.13692 0.60803 7.16470
H 30.04443 3.41265 7.41494
H 28.16282 1.42050 6.21869
H 28.08643 3.04243 5.49909
H 27.78412 2.88414 7.20297
H 32.31797 2.68554 7.46927
H 32.80296 0.96469 7.24868
H 30.86667 1.07551 9.33338
H 31.75626 2.54981 9.47464
H 30.15699 7.34692 5.60563
H 29.80780 8.48191 4.34506
H 31.47743 8.40796 4.95631
H 30.75798 9.70698 2.74139
H 32.02098 10.28193 1.50812
H 32.46212 9.39024 2.99260
H 32.01742 4.55767 -2.21692
H 30.76895 1.87398 -2.52666
H 30.10873 5.06634 -3.32475
H 29.21322 4.62846 -1.82730
H 29.45695 3.45813 -3.02168
H 32.77967 2.66075 -3.83426
H 32.91677 1.04583 -3.16479
H 34.14193 3.60336 -2.42985
H 34.95224 2.00176 -2.58182
H 33.98188 2.44447 -1.10644
H 32.50603 -3.25016 0.72030
H 31.99507 -2.19732 -0.63862
H 30.73655 -2.87098 0.47249
H 30.41475 -1.65100 5.44368
H 32.50668 -3.41381 8.77523
H 32.59861 -1.58221 9.17613
H 33.77740 -2.47887 8.18684

