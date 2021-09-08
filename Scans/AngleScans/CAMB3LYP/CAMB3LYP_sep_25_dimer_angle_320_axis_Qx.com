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
Mg 13.14598 1.18234 1.00561
C 11.33478 -0.47284 3.62458
C 13.56825 3.79753 3.15435
C 14.35041 2.84596 -1.41289
C 11.62773 -1.00863 -1.23945
N 12.49039 1.50244 3.09367
C 11.86150 0.77839 4.00723
C 11.86315 1.39679 5.43100
C 12.06799 2.93670 5.02409
C 12.78558 2.74355 3.68976
C 10.76229 3.66890 4.76261
C 13.05256 0.87167 6.30938
C 12.94070 1.19253 7.80146
H 13.43902 0.21320 8.76326
N 14.27008 2.81787 1.00559
C 14.33574 3.76014 1.96948
C 15.29201 4.76687 1.48349
C 15.66072 4.42583 0.08261
C 14.73865 3.33799 -0.17527
C 15.82485 5.86913 2.35214
C 16.66675 4.99418 -0.87189
O 16.65257 4.69994 -2.11931
C 17.61073 5.97742 -0.32086
N 12.85154 1.09539 -1.05062
C 13.53972 1.89728 -1.88206
C 13.36800 1.46858 -3.31981
C 12.61292 0.06005 -3.26087
C 12.38968 -0.04380 -1.73615
C 12.58234 2.42998 -4.21135
C 13.28911 -1.13209 -3.99179
C 14.57514 -1.50526 -3.29019
N 11.88299 -0.56466 1.15018
C 11.33948 -1.31483 0.08439
C 10.45444 -2.33753 0.67872
C 10.52709 -2.09474 2.05867
C 11.38023 -0.98535 2.30548
C 9.68255 -3.29086 -0.04341
C 10.06771 -2.51798 3.36195
O 9.39662 -3.45929 3.78301
C 10.47059 -1.40719 4.42507
C 11.09387 -2.09913 5.60628
O 12.17545 -2.70139 5.48762
O 10.42582 -1.75401 6.72206
C 11.01367 -2.40667 7.92218
H 13.56239 4.76003 3.66664
H 14.82564 3.31562 -2.27721
H 11.19493 -1.79157 -1.86584
H 10.93639 1.27720 5.97927
H 12.71251 3.45016 5.73283
H 9.91492 2.99967 4.77853
H 10.86911 4.12088 3.78519
H 10.64728 4.49773 5.46698
H 13.95837 1.44020 6.06234
H 13.18763 -0.18407 6.17593
H 11.93208 1.52339 8.07743
H 13.57479 2.05535 8.01935
H 15.24787 5.96444 3.28026
H 15.64354 6.81023 1.82174
H 16.89539 5.77295 2.54540
H 17.05704 6.81523 0.09623
H 18.30644 6.20101 -1.13195
H 18.15335 5.51369 0.51048
H 14.32914 1.33223 -3.80013
H 11.61554 0.11282 -3.72163
H 13.14843 2.75976 -5.10290
H 12.28297 3.28136 -3.61464
H 11.62653 2.05118 -4.56930
H 13.56295 -0.84579 -5.01022
H 12.65594 -2.01014 -4.06333
H 15.29702 -0.79462 -3.70370
H 14.85415 -2.53118 -3.52187
H 14.50528 -1.31495 -2.21290
H 9.79309 -4.23565 0.47201
H 10.00668 -3.36205 -1.07958
H 8.68839 -2.84483 0.04598
H 9.56148 -0.87002 4.69812
H 10.20550 -2.95543 8.41554
H 11.49574 -1.59440 8.49800
H 11.73712 -3.19692 7.75915

