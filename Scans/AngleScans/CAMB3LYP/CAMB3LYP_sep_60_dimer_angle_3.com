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
Mg 31.53890 2.84330 2.40253
C 31.34466 0.44703 5.06697
C 30.28278 5.14633 4.58262
C 31.27360 4.85844 -0.03391
C 31.67995 0.15866 0.18475
N 30.93208 2.69523 4.52434
C 30.96089 1.74893 5.45065
C 30.63269 2.24290 6.88517
C 29.77941 3.54790 6.50079
C 30.38810 3.84337 5.13153
C 28.29808 3.26487 6.31396
C 31.92155 2.61809 7.69763
C 31.70123 2.81220 9.19947
H 32.75800 2.39627 10.11746
N 31.34620 4.81842 2.38465
C 30.83660 5.59518 3.36354
C 30.89672 6.97433 2.85512
C 31.33069 6.93124 1.43227
C 31.31321 5.50181 1.19452
C 30.63672 8.17279 3.72088
C 31.68830 8.00040 0.44474
O 31.80690 7.74893 -0.80647
C 31.80384 9.36788 0.97214
N 31.27123 2.55914 0.36048
C 31.24162 3.60429 -0.48492
C 31.31716 3.14594 -1.92179
C 31.64873 1.58278 -1.85724
C 31.61799 1.38058 -0.32637
C 30.05568 3.36428 -2.75689
C 32.89716 1.09516 -2.64230
C 34.15366 1.64640 -2.00770
N 31.70408 0.69548 2.56883
C 31.72003 -0.24302 1.51382
C 31.72982 -1.58709 2.12687
C 31.69536 -1.33558 3.50674
C 31.64576 0.06586 3.73701
C 31.71836 -2.82331 1.42137
C 31.67892 -1.93722 4.82069
O 31.79158 -3.08373 5.25261
C 31.32344 -0.81316 5.88677
C 32.30141 -0.92570 7.02406
O 33.50967 -0.69260 6.84406
O 31.62238 -1.07597 8.17588
C 32.54892 -1.18098 9.33447
H 29.68388 5.88638 5.11414
H 31.29332 5.51183 -0.90922
H 31.82266 -0.72769 -0.43710
H 30.02759 1.56270 7.47250
H 29.97577 4.36528 7.18977
H 28.08161 2.20770 6.35481
H 28.04220 3.66618 5.34210
H 27.71092 3.83502 7.03957
H 32.23636 3.63261 7.42124
H 32.69733 1.89497 7.53670
H 30.73082 2.42043 9.52775
H 31.64141 3.88359 9.40580
H 30.17895 7.88731 4.67605
H 29.86750 8.76930 3.21856
H 31.52603 8.79048 3.86323
H 30.86186 9.65876 1.43097
H 32.15288 9.97545 0.13494
H 32.55664 9.37334 1.76830
H 32.11640 3.65327 -2.44804
H 30.83044 0.97526 -2.27095
H 30.23350 3.96865 -3.66638
H 29.30795 3.83161 -2.12996
H 29.55162 2.45454 -3.07861
H 32.87351 1.47652 -3.66609
H 32.97450 0.01472 -2.70281
H 34.22854 2.64910 -2.43907
H 35.01558 1.03736 -2.27267
H 34.02953 1.76191 -0.92469
H 32.43516 -3.46837 1.91187
H 31.96220 -2.68351 0.37029
H 30.67600 -3.12017 1.56497
H 30.29627 -0.98303 6.21179
H 32.30782 -2.11418 9.85282
H 32.42258 -0.24117 9.90417
H 33.60223 -1.32276 9.12262

