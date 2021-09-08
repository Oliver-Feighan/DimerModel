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
Mg 13.14655 1.18767 1.00308
C 13.57680 -0.97983 3.83078
C 12.49493 3.69969 3.21662
C 12.32122 3.00964 -1.45998
C 12.70253 -1.66889 -0.94336
N 13.07563 1.23075 3.21337
C 13.31701 0.36108 4.18278
C 13.35699 0.97918 5.60613
C 12.45403 2.27807 5.32968
C 12.71308 2.44190 3.83344
C 10.96821 2.03251 5.53310
C 14.81068 1.37321 6.04604
C 14.96750 1.69424 7.53403
H 16.21089 1.31664 8.20026
N 12.98315 3.16053 0.86281
C 12.73971 4.03023 1.86563
C 12.69294 5.36142 1.24122
C 12.76471 5.18975 -0.23529
C 12.66943 3.74736 -0.33790
C 12.66971 6.63352 2.03788
C 12.88469 6.16354 -1.36823
O 12.68980 5.80893 -2.58441
C 13.14507 7.56396 -1.00421
N 12.38311 0.75067 -0.88067
C 12.16215 1.72532 -1.78021
C 11.87715 1.15127 -3.14755
C 12.19240 -0.41228 -3.03131
C 12.53454 -0.49031 -1.52724
C 10.45282 1.34575 -3.66701
C 13.20364 -1.00405 -4.05105
C 14.58490 -0.44772 -3.79027
N 13.31715 -0.94442 1.30830
C 13.06112 -1.96422 0.36568
C 13.20177 -3.25452 1.07120
C 13.50972 -2.89250 2.39128
C 13.53774 -1.47627 2.50535
C 13.00053 -4.54203 0.49854
C 13.80700 -3.38628 3.71661
O 14.00581 -4.49782 4.20516
C 13.73919 -2.16899 4.73641
C 14.96415 -2.22410 5.60759
O 16.09472 -2.04813 5.12031
O 14.58572 -2.25822 6.89833
C 15.76614 -2.30234 7.80193
H 12.05486 4.50016 3.81199
H 12.13525 3.58991 -2.36672
H 12.67616 -2.60650 -1.50270
H 12.90458 0.36947 6.37901
H 12.82458 3.14060 5.87749
H 10.75346 0.99012 5.71599
H 10.48786 2.36350 4.62164
H 10.58465 2.67883 6.32777
H 15.06250 2.35088 5.61535
H 15.51320 0.61311 5.76437
H 14.10158 1.36381 8.12049
H 14.97511 2.78012 7.65583
H 12.45574 6.44132 3.09650
H 11.80939 7.21427 1.68796
H 13.57540 7.22942 1.90644
H 12.34828 7.92304 -0.35719
H 13.28705 8.09023 -1.95013
H 14.06988 7.60695 -0.41813
H 12.53029 1.58692 -3.89369
H 11.28925 -1.02220 -3.18020
H 10.41106 1.86898 -4.64095
H 9.88798 1.88731 -2.91988
H 9.87258 0.43122 -3.77731
H 12.93544 -0.70513 -5.06735
H 13.24860 -2.08790 -4.03542
H 14.56600 0.51410 -4.31146
H 15.34706 -1.10549 -4.20320
H 14.73132 -0.24180 -2.72359
H 13.80645 -5.17027 0.85393
H 12.98157 -4.49520 -0.58830
H 12.02097 -4.79017 0.91586
H 12.82056 -2.27664 5.31443
H 15.64616 -3.18222 8.44154
H 15.79636 -1.31620 8.30230
H 16.73343 -2.49701 7.35371

