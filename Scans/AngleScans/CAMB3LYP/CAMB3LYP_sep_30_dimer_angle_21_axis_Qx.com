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
Mg 15.77055 1.42656 1.20421
C 16.46993 -0.61274 4.07322
C 13.76228 3.30941 3.21825
C 14.78317 2.96699 -1.38800
C 16.88050 -1.21930 -0.77396
N 15.26295 1.26901 3.35125
C 15.63716 0.49019 4.35511
C 15.15309 0.96401 5.75180
C 13.88110 1.82635 5.28565
C 14.33761 2.19228 3.87485
C 12.60564 1.00729 5.17756
C 16.21660 1.85754 6.48157
C 15.94269 2.09915 7.96762
H 17.07967 2.18512 8.87996
N 14.86882 3.18489 1.02063
C 14.11184 3.81168 1.94551
C 13.66268 5.06345 1.31662
C 14.08033 5.04700 -0.11167
C 14.58649 3.69366 -0.22274
C 12.98363 6.16113 2.08304
C 14.02074 6.07453 -1.20110
O 14.22130 5.76644 -2.42897
C 13.62890 7.43353 -0.79977
N 15.62243 0.87285 -0.79424
C 15.21144 1.75063 -1.72618
C 15.44732 1.21762 -3.11925
C 16.32770 -0.10403 -2.93021
C 16.37523 -0.15791 -1.38741
C 14.19219 0.88189 -3.92444
C 17.66686 -0.17568 -3.71392
C 18.63562 0.85240 -3.17544
N 16.71003 -0.48779 1.55164
C 17.06661 -1.45131 0.58284
C 17.56814 -2.63523 1.31032
C 17.44605 -2.28433 2.66326
C 16.88768 -0.98204 2.77173
C 18.00859 -3.85157 0.71647
C 17.65266 -2.72362 4.02452
O 18.17744 -3.70416 4.55071
C 16.91223 -1.71070 5.00030
C 17.86523 -1.35156 6.10724
O 18.90428 -0.71282 5.86432
O 17.28984 -1.62909 7.29142
C 18.19224 -1.27967 8.42075
H 12.93497 3.82740 3.70417
H 14.56132 3.49689 -2.31723
H 17.33660 -2.04742 -1.32065
H 14.83951 0.16901 6.41769
H 13.76591 2.72032 5.89303
H 12.79085 -0.04722 5.31843
H 12.21933 1.19418 4.18414
H 11.85167 1.39047 5.87114
H 16.13817 2.88606 6.10642
H 17.20287 1.45465 6.35675
H 15.18326 1.41396 8.36350
H 15.49548 3.08963 8.08150
H 12.66334 5.82026 3.07533
H 12.04882 6.38642 1.55842
H 13.58562 7.07061 2.13884
H 12.64646 7.40358 -0.33463
H 13.73037 8.04444 -1.69887
H 14.32869 7.78783 -0.03466
H 16.00490 1.92888 -3.71618
H 15.78770 -1.00389 -3.25959
H 14.13537 1.42068 -4.88909
H 13.32624 1.10202 -3.31435
H 14.05531 -0.17494 -4.14709
H 17.50391 0.07234 -4.76563
H 18.13390 -1.15439 -3.68246
H 18.33798 1.76799 -3.69524
H 19.66018 0.57672 -3.41690
H 18.47938 1.01660 -2.10293
H 18.91234 -4.14204 1.23546
H 18.18291 -3.73261 -0.35088
H 17.14721 -4.49245 0.92279
H 16.01881 -2.21124 5.37558
H 18.30988 -2.18315 9.02703
H 17.73176 -0.40099 8.91026
H 19.22414 -1.04771 8.18445

