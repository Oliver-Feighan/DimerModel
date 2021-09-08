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
Mg 7.89124 0.70445 0.59446
C 7.81996 2.01041 -2.74747
C 7.38921 -2.38473 -0.76182
C 7.49947 -0.41475 3.53747
C 7.24530 3.93779 1.72317
N 7.65210 0.02816 -1.49784
C 7.70451 0.60590 -2.68856
C 7.70971 -0.38411 -3.88414
C 7.00256 -1.63362 -3.16475
C 7.39749 -1.33897 -1.71911
C 5.48624 -1.60440 -3.26108
C 9.16221 -0.74095 -4.35839
C 9.24035 -1.45446 -5.70988
H 10.36965 -1.16059 -6.58784
N 7.98973 -1.15325 1.28569
C 7.77960 -2.29174 0.59217
C 7.95064 -3.38839 1.55771
C 8.11649 -2.79793 2.91365
C 7.84834 -1.40269 2.62818
C 8.02534 -4.83022 1.14665
C 8.44746 -3.39270 4.24891
O 8.30598 -2.72818 5.33575
C 8.85294 -4.80586 4.25478
N 7.22985 1.58096 2.35986
C 7.20536 0.88528 3.51023
C 6.95901 1.79311 4.69156
C 7.06424 3.28234 4.11825
C 7.27400 2.96317 2.62175
C 5.61609 1.61423 5.39934
C 8.07003 4.23663 4.81854
C 9.48553 3.76876 4.56786
N 7.76677 2.66783 -0.29868
C 7.45929 3.88290 0.35182
C 7.37984 4.92679 -0.69049
C 7.62586 4.23687 -1.88714
C 7.82330 2.85666 -1.61231
C 7.06356 6.29743 -0.47271
C 7.75290 4.36052 -3.32146
O 7.77066 5.30137 -4.11407
C 7.75898 2.90312 -3.95575
C 8.89457 2.82938 -4.93938
O 10.07322 2.91202 -4.55135
O 8.41418 2.45761 -6.13995
C 9.50455 2.35976 -7.14659
H 7.00807 -3.36088 -1.06309
H 7.46037 -0.72883 4.58308
H 7.14504 4.98829 2.00437
H 7.12439 -0.06713 -4.73895
H 7.43459 -2.57516 -3.49371
H 5.12790 -0.68359 -3.69677
H 5.12496 -1.70856 -2.24635
H 5.12580 -2.48458 -3.80117
H 9.56860 -1.52595 -3.70772
H 9.78320 0.13353 -4.37500
H 8.29608 -1.39178 -6.26429
H 9.37530 -2.52386 -5.53042
H 7.70594 -4.96899 0.10639
H 7.27552 -5.36979 1.73528
H 9.00653 -5.27145 1.33440
H 8.05928 -5.41109 3.82320
H 9.13442 -5.02483 5.28662
H 9.72661 -4.92185 3.60373
H 7.71875 1.65439 5.45094
H 6.10593 3.81686 4.19424
H 5.71771 1.38792 6.47755
H 5.06692 0.82951 4.89620
H 4.93561 2.46025 5.31935
H 7.92275 4.21385 5.90111
H 7.97646 5.27023 4.50248
H 9.62936 2.99779 5.33075
H 10.18877 4.58885 4.69799
H 9.57211 3.28372 3.58862
H 7.75316 6.87534 -1.07352
H 7.13650 6.55978 0.58067
H 6.03066 6.31846 -0.83040
H 6.79138 2.75054 -4.43538
H 9.22438 3.00518 -7.98482
H 9.61940 1.28013 -7.35862
H 10.47188 2.76858 -6.87901

