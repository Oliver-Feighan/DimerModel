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
Mg 10.50077 0.94884 0.79640
C 9.46503 -1.86956 -1.16901
C 9.75194 2.90956 -1.89235
C 11.82311 3.38945 2.33339
C 12.20610 -1.28856 2.85349
N 9.71604 0.46581 -1.21435
C 9.30114 -0.63277 -1.82698
C 8.60183 -0.38388 -3.19027
C 9.26631 1.03414 -3.54604
C 9.56872 1.52415 -2.13144
C 10.58004 0.91094 -4.29989
C 7.04323 -0.26892 -3.05037
C 6.27243 -0.33487 -4.37070
H 4.97021 -0.99555 -4.38418
N 10.31213 2.90309 0.50673
C 9.97910 3.52691 -0.64271
C 9.99585 4.96619 -0.33880
C 10.52624 5.14839 1.03968
C 10.93619 3.79253 1.34574
C 9.46275 6.00363 -1.28375
C 10.65624 6.34722 1.92969
O 11.36473 6.31948 2.99744
C 10.00601 7.58158 1.46597
N 11.99307 1.06467 2.23932
C 12.34222 2.24542 2.77953
C 13.23294 2.05233 3.98358
C 13.21268 0.48198 4.28556
C 12.34712 0.00704 3.09786
C 14.67903 2.52086 3.82327
C 12.80347 0.04623 5.71921
C 11.34638 0.36877 5.96052
N 10.65164 -1.19946 0.96980
C 11.44275 -1.92789 1.88502
C 11.30225 -3.35753 1.54036
C 10.45064 -3.36005 0.42535
C 10.10226 -2.02482 0.08581
C 11.95606 -4.43923 2.19503
C 9.77568 -4.18797 -0.54813
O 9.62981 -5.39934 -0.70640
C 9.20884 -3.25224 -1.70117
C 7.78167 -3.64927 -1.96207
O 6.91060 -3.48594 -1.08961
O 7.63964 -3.95163 -3.26544
C 6.23822 -4.34099 -3.57567
H 9.76515 3.58628 -2.74720
H 12.22127 4.18778 2.96394
H 12.62774 -2.05517 3.50710
H 8.83720 -1.10850 -3.96053
H 8.55289 1.69345 -4.03358
H 10.91040 -0.11486 -4.36868
H 11.29522 1.50223 -3.74318
H 10.49687 1.38481 -5.28217
H 6.78636 0.75009 -2.73354
H 6.66564 -1.01099 -2.37406
H 6.89717 -0.70001 -5.19491
H 6.00481 0.68215 -4.66745
H 9.29159 5.58673 -2.28390
H 10.25846 6.74213 -1.42984
H 8.57667 6.51085 -0.89613
H 10.40907 7.85924 0.49493
H 10.13340 8.30139 2.27685
H 8.93906 7.38067 1.31807
H 12.84040 2.58068 4.84372
H 14.20580 0.02849 4.15180
H 14.98369 3.26460 4.58361
H 14.79841 2.92500 2.82680
H 15.42586 1.72899 3.84241
H 13.37316 0.60863 6.46306
H 12.97096 -1.00768 5.91465
H 11.37166 1.42392 6.24882
H 10.94460 -0.24952 6.76072
H 10.76788 0.29442 5.03237
H 11.21681 -5.21888 2.32205
H 12.37714 -4.13023 3.14949
H 12.73153 -4.67184 1.46031
H 9.84041 -3.39475 -2.57889
H 6.27723 -5.33044 -4.04166
H 5.82594 -3.51381 -4.18356
H 5.57304 -4.52663 -2.74049

