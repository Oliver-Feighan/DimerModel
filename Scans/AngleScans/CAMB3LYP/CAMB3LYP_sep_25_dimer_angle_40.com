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
Mg 13.14526 1.18262 1.00627
C 13.19740 -2.32822 1.74812
C 11.97667 1.68201 4.17164
C 12.66066 4.22351 0.24101
C 13.20967 0.32833 -2.37186
N 12.68756 -0.22109 2.65319
C 12.80487 -1.52703 2.84057
C 12.56097 -1.99165 4.30146
C 11.64817 -0.76583 4.79444
C 12.15457 0.31478 3.84129
C 10.16593 -0.97630 4.53368
C 13.89102 -2.08634 5.12851
C 13.76749 -2.82576 6.46259
H 14.89492 -3.63514 6.91678
N 12.89766 2.77476 2.16487
C 12.43428 2.79474 3.43222
C 12.42203 4.20992 3.83384
C 12.75984 5.03803 2.64437
C 12.76556 4.02332 1.60975
C 12.18837 4.65208 5.24922
C 13.02041 6.50303 2.46527
O 13.06100 7.04261 1.30336
C 13.13377 7.30365 3.69313
N 12.74805 2.13650 -0.79780
C 12.63269 3.47462 -0.86158
C 12.62347 3.95369 -2.29364
C 13.00141 2.67513 -3.17683
C 13.07980 1.61066 -2.06052
C 11.30291 4.54225 -2.78973
C 14.20683 2.82037 -4.14562
C 15.48775 2.97022 -3.35685
N 13.38017 -0.63776 -0.13351
C 13.35051 -0.77413 -1.53872
C 13.43860 -2.21813 -1.83794
C 13.49061 -2.82805 -0.57549
C 13.41825 -1.83570 0.43915
C 13.41340 -2.80176 -3.13596
C 13.57955 -4.08647 0.12962
O 13.75257 -5.25828 -0.20303
C 13.26624 -3.82794 1.66631
C 14.32158 -4.52627 2.47940
O 15.50805 -4.15751 2.42539
O 13.72634 -5.36635 3.34562
C 14.73162 -6.07430 4.18230
H 11.39499 1.92968 5.06001
H 12.60323 5.26617 -0.07998
H 13.33428 -0.01259 -3.40187
H 12.01585 -2.92306 4.39691
H 11.86816 -0.49896 5.82485
H 9.98180 -1.86668 3.95104
H 9.83401 -0.09739 3.99665
H 9.61374 -0.97901 5.47782
H 14.15847 -1.08589 5.49226
H 14.67370 -2.52711 4.54219
H 12.83261 -3.39491 6.53338
H 13.69238 -2.08627 7.26355
H 11.80422 3.83217 5.86857
H 11.37087 5.38078 5.22575
H 13.06798 5.12197 5.69437
H 12.21733 7.21035 4.27098
H 13.40861 8.30738 3.36319
H 13.93823 6.88699 4.30945
H 13.37108 4.72177 -2.44931
H 12.17401 2.37728 -3.83753
H 11.40219 5.57526 -3.17329
H 10.58673 4.50460 -1.97966
H 10.80341 3.96600 -3.56667
H 14.10356 3.72810 -4.74522
H 14.30952 1.98899 -4.83501
H 15.50575 4.03731 -3.11584
H 16.34617 2.68796 -3.96306
H 15.43399 2.41942 -2.41061
H 14.17914 -3.56594 -3.14865
H 13.58170 -2.05631 -3.91040
H 12.39166 -3.19029 -3.15531
H 12.26844 -4.21966 1.86762
H 14.55180 -7.14673 4.05883
H 14.61833 -5.65859 5.20118
H 15.77169 -5.99890 3.88733

