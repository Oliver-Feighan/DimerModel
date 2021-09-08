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
Mg 15.77298 1.42744 1.19599
C 15.16353 3.43421 4.10812
C 14.79054 4.10665 -0.67249
C 15.86713 -0.44778 -1.36177
C 15.59926 -1.37455 3.26103
N 15.08874 3.45372 1.76013
C 14.95620 4.12347 2.89512
C 14.65867 5.63664 2.71806
C 14.00994 5.58397 1.24998
C 14.69454 4.32354 0.72531
C 12.50968 5.34232 1.26329
C 15.96294 6.50868 2.73875
C 15.72783 8.01407 2.88125
H 16.68763 8.79842 3.65328
N 15.86166 1.88886 -0.73217
C 15.43805 3.02933 -1.31626
C 15.70731 2.86926 -2.75364
C 16.17282 1.47588 -2.99135
C 15.96076 0.89956 -1.67873
C 15.59357 4.00059 -3.73373
C 16.70689 0.77433 -4.20329
O 16.82581 -0.50118 -4.24522
C 16.99865 1.61618 -5.37276
N 15.52849 -0.62121 0.93988
C 15.67179 -1.18855 -0.27070
C 15.72440 -2.69401 -0.16481
C 15.82992 -3.01126 1.39909
C 15.72560 -1.57461 1.95635
C 14.53169 -3.44461 -0.75677
C 17.01971 -3.89775 1.85867
C 18.32209 -3.15553 1.66297
N 15.62828 1.06841 3.32160
C 15.54283 -0.18233 3.97153
C 15.34460 0.08736 5.41042
C 15.30512 1.48726 5.49533
C 15.44707 2.04992 4.19811
C 15.17964 -0.89613 6.42613
C 15.16499 2.61675 6.38601
O 15.10181 2.75792 7.60659
C 14.94033 3.92497 5.51167
C 15.85849 4.99592 6.03376
O 17.09234 4.87106 5.94093
O 15.13121 6.08089 6.35669
C 15.99891 7.17460 6.86924
H 14.28662 4.80465 -1.34164
H 16.00483 -1.13916 -2.19634
H 15.63377 -2.19309 3.98316
H 13.94647 6.04531 3.42484
H 14.29922 6.44911 0.65905
H 12.14501 5.12744 2.25677
H 12.34192 4.49779 0.60780
H 11.98765 6.18748 0.80555
H 16.42587 6.48417 1.74384
H 16.63329 6.17341 3.50611
H 14.70254 8.24312 3.19628
H 15.81378 8.47373 1.89365
H 15.07198 4.86083 -3.29632
H 14.93149 3.66166 -4.53787
H 16.55656 4.28325 -4.16433
H 16.09399 2.13688 -5.67762
H 17.45462 0.94891 -6.10672
H 17.72075 2.38592 -5.07832
H 16.60250 -3.08628 -0.66299
H 14.94667 -3.55516 1.76501
H 14.81992 -4.18188 -1.52970
H 13.83918 -2.71919 -1.16265
H 13.91400 -3.97406 -0.03334
H 17.08048 -4.79876 1.24333
H 16.94550 -4.21821 2.89249
H 18.55052 -3.33197 0.60762
H 19.09686 -3.56464 2.30827
H 18.18347 -2.07608 1.79423
H 15.78334 -0.58012 7.26647
H 15.47181 -1.88350 6.07469
H 14.10201 -0.82291 6.59533
H 13.89029 4.20486 5.60472
H 15.61323 7.45316 7.85480
H 15.98969 7.95538 6.08562
H 17.02747 6.92901 7.10632

