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
Mg 10.51626 0.95240 0.80519
C 10.76394 -1.22243 3.64910
C 8.85737 3.15024 2.81819
C 9.88672 2.68747 -1.77563
C 11.14915 -1.82655 -1.20046
N 9.94580 0.86259 2.94031
C 10.14739 0.01319 3.93642
C 9.73329 0.54967 5.33293
C 8.65506 1.64211 4.86083
C 9.19989 1.93573 3.46454
C 7.25101 1.07884 4.71698
C 10.93088 1.21709 6.09605
C 10.67784 1.48476 7.58134
H 11.79252 1.34269 8.51395
N 9.96411 2.85124 0.63756
C 9.32000 3.59589 1.56054
C 9.12610 4.91853 0.94638
C 9.56138 4.84417 -0.47478
C 9.80684 3.42155 -0.60120
C 8.64996 6.11326 1.72056
C 9.71710 5.87994 -1.54677
O 9.88054 5.55709 -2.77636
C 9.57925 7.28253 -1.12816
N 10.30643 0.46475 -1.20492
C 10.08585 1.41723 -2.12799
C 10.24506 0.86920 -3.52610
C 10.85803 -0.59685 -3.34567
C 10.86412 -0.68056 -1.80348
C 8.96540 0.78666 -4.35798
C 12.17521 -0.90765 -4.10798
C 13.30876 -0.08759 -3.53532
N 11.07298 -1.10908 1.13449
C 11.26159 -2.10859 0.15497
C 11.51771 -3.37584 0.86975
C 11.43689 -3.02746 2.22643
C 11.13065 -1.64513 2.34849
C 11.73386 -4.64472 0.26201
C 11.53052 -3.51699 3.58301
O 11.85158 -4.58601 4.10054
C 10.97407 -2.39695 4.56396
C 11.95546 -2.23894 5.69286
O 13.10045 -1.80337 5.47867
O 11.31494 -2.42020 6.86216
C 12.24438 -2.26256 8.01243
H 8.13245 3.80747 3.29925
H 9.78659 3.26272 -2.69896
H 11.45255 -2.71776 -1.75396
H 9.26308 -0.18163 5.97936
H 8.69760 2.53311 5.48186
H 7.23234 0.00642 4.84231
H 6.92631 1.34902 3.72068
H 6.56873 1.58695 5.40448
H 11.05417 2.24724 5.73783
H 11.82635 0.63791 5.98068
H 9.79570 0.94888 7.95231
H 10.42216 2.53990 7.70515
H 8.25188 5.82459 2.70118
H 7.78454 6.51755 1.18441
H 9.41062 6.89258 1.80245
H 8.59962 7.43107 -0.68014
H 9.81123 7.87618 -2.01453
H 10.31780 7.48819 -0.34531
H 10.93782 1.47144 -4.10097
H 10.16545 -1.37449 -3.69986
H 9.02970 1.34015 -5.31380
H 8.14422 1.15686 -3.75869
H 8.63716 -0.22239 -4.60149
H 12.08248 -0.61856 -5.15774
H 12.44969 -1.95702 -4.08597
H 13.19845 0.87489 -4.04381
H 14.26801 -0.54737 -3.76442
H 13.16493 0.08782 -2.46285
H 12.55667 -5.10711 0.79086
H 11.94845 -4.54550 -0.80000
H 10.76367 -5.11525 0.44255
H 9.99539 -2.72602 4.91534
H 12.17847 -3.18055 8.60457
H 11.94728 -1.32005 8.50957
H 13.30593 -2.22523 7.79757

