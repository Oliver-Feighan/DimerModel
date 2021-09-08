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
Mg 7.88895 0.70547 0.60662
C 5.35786 0.05697 3.06671
C 9.39338 2.72792 2.90466
C 9.88183 1.64878 -1.67505
C 5.63549 -0.41339 -1.80497
N 7.34966 1.24226 2.68346
C 6.39871 0.88629 3.53394
C 6.61599 1.38563 4.98756
C 7.55951 2.64581 4.67043
C 8.17086 2.16845 3.35473
C 6.78468 3.92497 4.40074
C 7.35332 0.32609 5.87960
C 7.32422 0.61796 7.38151
H 7.23368 -0.50746 8.30764
N 9.65924 1.59340 0.73403
C 10.11475 2.35844 1.74819
C 11.46356 2.79124 1.35125
C 11.70238 2.35497 -0.05137
C 10.38717 1.85608 -0.39982
C 12.41017 3.47325 2.29581
C 12.91067 2.39363 -0.93724
O 12.82827 2.17892 -2.19817
C 14.17789 2.78198 -0.30098
N 7.70758 0.83085 -1.46068
C 8.74314 1.22391 -2.22294
C 8.46916 0.97316 -3.68665
C 7.12752 0.10384 -3.72984
C 6.79476 0.07695 -2.22188
C 8.29570 2.21972 -4.55405
C 7.18733 -1.24502 -4.49774
C 8.09245 -2.21319 -3.77056
N 5.93524 -0.21705 0.61251
C 5.16002 -0.58015 -0.51056
C 3.85962 -1.06446 -0.00376
C 3.96079 -0.92685 1.38883
C 5.22672 -0.37547 1.72487
C 2.76714 -1.50486 -0.80292
C 3.28049 -1.11272 2.65032
O 2.21635 -1.62423 2.99629
C 4.10638 -0.36604 3.78468
C 4.25084 -1.30669 4.94951
O 4.91395 -2.35292 4.83857
O 3.76849 -0.71401 6.05695
C 3.89977 -1.60371 7.24144
H 9.82199 3.55820 3.46662
H 10.57314 1.85498 -2.49545
H 4.91612 -0.87131 -2.48724
H 5.71631 1.71303 5.49467
H 8.32968 2.76324 5.42832
H 5.72052 3.74840 4.35039
H 7.15171 4.29698 3.45314
H 7.04247 4.68552 5.14333
H 8.43341 0.39298 5.69610
H 6.97113 -0.65913 5.69548
H 6.58530 1.38691 7.63774
H 8.28156 1.06026 7.66768
H 11.89806 3.80875 3.20608
H 12.73486 4.39962 1.80972
H 13.28925 2.86622 2.52253
H 14.07275 3.77063 0.13944
H 14.94059 2.66478 -1.07326
H 14.38169 2.09006 0.52384
H 9.27188 0.40323 -4.13824
H 6.30692 0.64524 -4.22313
H 9.00074 2.26040 -5.40573
H 8.40925 3.09219 -3.92441
H 7.29789 2.36044 -4.96608
H 7.62309 -1.09751 -5.48895
H 6.21508 -1.70531 -4.63863
H 9.08958 -1.92830 -4.11934
H 7.85615 -3.23905 -4.04589
H 8.06108 -2.04392 -2.68796
H 2.37964 -2.39983 -0.33454
H 3.07589 -1.69396 -1.82893
H 2.10701 -0.63675 -0.72628
H 3.55395 0.53543 4.05252
H 2.90046 -1.70756 7.67534
H 4.67926 -1.14248 7.87658
H 4.16171 -2.64009 7.06324

