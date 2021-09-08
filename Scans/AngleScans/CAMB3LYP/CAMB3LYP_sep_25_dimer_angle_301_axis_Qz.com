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
Mg 13.14509 1.19118 0.99831
C 12.57598 2.39514 4.33082
C 12.07350 4.22943 -0.12209
C 13.20582 0.01505 -1.94795
C 13.06789 -2.04013 2.30153
N 12.44255 2.99754 2.06403
C 12.32749 3.36025 3.33271
C 12.00133 4.86296 3.54499
C 11.31724 5.16368 2.12357
C 12.00880 4.08899 1.28723
C 9.82186 4.89370 2.10622
C 13.29151 5.73054 3.75657
C 13.03580 7.14742 4.27504
H 14.00174 7.73536 5.19904
N 13.17838 2.12037 -0.75499
C 12.72207 3.36091 -1.02727
C 12.95806 3.57000 -2.46420
C 13.43989 2.29031 -3.05151
C 13.26984 1.40058 -1.92034
C 12.80180 4.90713 -3.12835
C 13.95490 1.92473 -4.41067
O 14.09327 0.70282 -4.77207
C 14.20391 3.03761 -5.33856
N 12.92733 -0.73384 0.24384
C 13.04961 -0.97833 -1.07267
C 13.12907 -2.46116 -1.34703
C 13.27852 -3.15577 0.08572
C 13.16497 -1.90594 0.98595
C 11.93425 -3.06637 -2.08358
C 14.49349 -4.10267 0.28548
C 15.77848 -3.30689 0.25520
N 13.05907 0.31070 2.96945
C 13.00995 -1.06411 3.28815
C 12.84323 -1.16589 4.75248
C 12.78334 0.16748 5.18494
C 12.88392 1.03862 4.06669
C 12.71943 -2.37483 5.49361
C 12.64723 1.03604 6.33202
O 12.61214 0.86716 7.55021
C 12.37987 2.51566 5.81663
C 13.29336 3.44246 6.57102
O 14.52638 3.37159 6.42524
O 12.55696 4.39660 7.16910
C 13.41941 5.34672 7.92094
H 11.54194 5.06101 -0.58556
H 13.33384 -0.44336 -2.93129
H 13.13352 -3.01185 2.79562
H 11.30043 5.06697 4.34555
H 11.57778 6.15488 1.76164
H 9.48551 4.43012 3.02170
H 9.65149 4.23572 1.26415
H 9.27508 5.81468 1.88451
H 13.72990 5.96483 2.77797
H 13.98603 5.22930 4.40234
H 12.01510 7.26834 4.65778
H 13.08976 7.84046 3.43195
H 12.27746 5.61965 -2.47967
H 12.12548 4.76495 -3.97820
H 13.74911 5.30910 -3.49398
H 13.28369 3.59801 -5.48558
H 14.65218 2.58440 -6.22486
H 14.92062 3.72523 -4.87584
H 14.00072 -2.69767 -1.94486
H 12.41352 -3.79281 0.32190
H 12.21502 -3.58130 -3.02172
H 11.22028 -2.27796 -2.28155
H 11.34335 -3.77273 -1.50300
H 14.55344 -4.82036 -0.53639
H 14.45017 -4.67219 1.20785
H 15.98342 -3.20975 -0.81515
H 16.57554 -3.84692 0.76226
H 15.62580 -2.29750 0.65453
H 13.33866 -2.26513 6.37395
H 13.01863 -3.23686 4.90103
H 11.64530 -2.36959 5.69731
H 11.32809 2.74057 5.99764
H 13.05390 5.36241 8.95237
H 13.37814 6.29775 7.35740
H 14.45736 5.07227 8.06853

