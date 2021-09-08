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
Mg 31.53840 2.84372 2.40268
C 31.15741 0.60456 5.18118
C 30.13973 5.26583 4.35503
C 31.44460 4.71453 -0.15832
C 31.82653 0.03605 0.35825
N 30.78710 2.81685 4.48284
C 30.75054 1.92557 5.46193
C 30.32554 2.50067 6.83968
C 29.50300 3.77951 6.32307
C 30.20478 3.99690 4.98407
C 28.03749 3.48297 6.05164
C 31.55626 2.92485 7.71571
C 31.23371 3.20471 9.18520
H 32.22422 2.84466 10.19600
N 31.35093 4.81411 2.25780
C 30.77674 5.64488 3.15300
C 30.87409 6.99257 2.57126
C 31.40464 6.86849 1.18640
C 31.40095 5.42771 1.03081
C 30.55744 8.23836 3.34662
C 31.83112 7.87974 0.16577
O 32.03487 7.55685 -1.05777
C 31.91264 9.27558 0.62015
N 31.41102 2.44178 0.36686
C 31.44140 3.43641 -0.53742
C 31.61457 2.89618 -1.93690
C 31.93812 1.34006 -1.75976
C 31.80200 1.22635 -0.22550
C 30.41379 3.06330 -2.86781
C 33.23662 0.81071 -2.42786
C 34.44758 1.40033 -1.74135
N 31.68791 0.70942 2.70344
C 31.77456 -0.28829 1.70787
C 31.73982 -1.59476 2.39659
C 31.61117 -1.26423 3.75405
C 31.54841 0.14803 3.89921
C 31.77460 -2.86961 1.76440
C 31.50350 -1.78919 5.09626
O 31.58418 -2.90866 5.60026
C 31.07770 -0.60634 6.06886
C 31.97509 -0.65104 7.27518
O 33.19326 -0.42608 7.16540
O 31.21832 -0.73615 8.38443
C 32.06295 -0.77221 9.60799
H 29.50709 6.03398 4.80066
H 31.52554 5.31644 -1.06641
H 32.01000 -0.88435 -0.20022
H 29.68033 1.85414 7.42233
H 29.65307 4.63565 6.97564
H 27.81682 2.42945 6.13839
H 27.84964 3.82706 5.04299
H 27.40293 4.09272 6.70116
H 31.89112 3.92243 7.40356
H 32.33995 2.19535 7.65035
H 30.24235 2.83040 9.46813
H 31.16179 4.28607 9.32484
H 30.03466 8.00741 4.28298
H 29.82558 8.80327 2.75914
H 31.43598 8.86515 3.51383
H 30.94191 9.59039 0.99572
H 32.31942 9.83466 -0.22474
H 32.60901 9.32855 1.46447
H 32.44896 3.37407 -2.43542
H 31.14907 0.70794 -2.19295
H 30.65471 3.61465 -3.79627
H 29.62563 3.56436 -2.32166
H 29.93137 2.13545 -3.17039
H 33.28398 1.13239 -3.47114
H 33.31599 -0.27125 -2.42055
H 34.55369 2.37666 -2.22364
H 35.32455 0.77890 -1.91099
H 34.24960 1.57779 -0.67788
H 32.45488 -3.48378 2.33931
H 32.09027 -2.79008 0.72622
H 30.72431 -3.15996 1.85299
H 30.03033 -0.75942 6.33187
H 31.78516 -1.67451 10.16150
H 31.89949 0.19859 10.11256
H 33.12807 -0.92368 9.47747

