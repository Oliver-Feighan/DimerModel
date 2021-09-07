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
Mg 8.94146 0.79975 0.68492
C 9.13995 -2.46552 -0.79079
C 7.97627 -0.72954 3.57702
C 8.32298 3.67277 1.88103
C 8.80078 2.09920 -2.54577
N 8.63483 -1.32880 1.20214
C 8.80033 -2.48529 0.57804
C 8.66859 -3.73359 1.49124
C 7.75791 -3.09081 2.64726
C 8.16887 -1.62753 2.49693
C 6.26777 -3.19935 2.36973
C 10.05361 -4.21341 2.05117
C 10.04088 -5.60145 2.69532
H 11.21821 -6.45140 2.54060
N 8.72916 1.38826 2.56831
C 8.35214 0.63105 3.61986
C 8.32824 1.53613 4.77936
C 8.56202 2.92406 4.29592
C 8.52563 2.71358 2.86258
C 8.17875 1.04681 6.19061
C 8.76958 4.22745 5.00606
O 8.71672 5.34796 4.38600
C 8.94368 4.15905 6.46432
N 8.39718 2.60466 -0.19165
C 8.24110 3.71461 0.55104
C 8.12195 4.94193 -0.32067
C 8.47426 4.45269 -1.80191
C 8.65711 2.94284 -1.53295
C 6.75522 5.62635 -0.32262
C 9.60699 5.21382 -2.54373
C 10.93371 4.95075 -1.86846
N 9.14864 0.01523 -1.31750
C 9.02782 0.72895 -2.52997
C 9.13506 -0.25484 -3.62682
C 9.28901 -1.48544 -2.97075
C 9.25821 -1.28539 -1.56415
C 9.03819 0.03533 -5.01698
C 9.45991 -2.90951 -3.14823
O 9.64213 -3.64879 -4.11467
C 9.24422 -3.62320 -1.74448
C 10.37088 -4.59889 -1.54148
O 11.54041 -4.19605 -1.41299
O 9.85883 -5.82252 -1.31589
C 10.93742 -6.82294 -1.09762
H 7.44939 -1.08773 4.46189
H 8.21540 4.69904 2.23972
H 8.86478 2.43702 -3.58233
H 8.15693 -4.57495 1.03948
H 8.03965 -3.46724 3.62720
H 6.06917 -3.58684 1.38152
H 5.87632 -2.19529 2.46813
H 5.78095 -3.79049 3.15068
H 10.31748 -3.60334 2.92464
H 10.80665 -4.17601 1.28819
H 9.12889 -6.16023 2.45245
H 9.99982 -5.48005 3.78049
H 7.85999 -0.00241 6.22093
H 7.34185 1.59784 6.63314
H 9.07322 1.21896 6.79308
H 8.07131 3.68743 6.91043
H 9.16795 5.17978 6.78005
H 9.79910 3.51066 6.68450
H 8.83597 5.69946 -0.02163
H 7.61257 4.54923 -2.47867
H 6.79996 6.69132 -0.02610
H 6.09674 5.07562 0.33590
H 6.22034 5.58692 -1.27001
H 9.43855 6.29207 -2.48756
H 9.68561 4.95444 -3.59422
H 10.93868 5.67122 -1.04506
H 11.75660 5.13242 -2.55685
H 10.95917 3.94707 -1.42841
H 9.82201 -0.52629 -5.50761
H 9.13323 1.10221 -5.20800
H 8.02853 -0.33016 -5.22261
H 8.27344 -4.11941 -1.77498
H 10.77916 -7.62682 -1.82318
H 10.88191 -7.09317 -0.02634
H 11.95271 -6.52430 -1.33095

