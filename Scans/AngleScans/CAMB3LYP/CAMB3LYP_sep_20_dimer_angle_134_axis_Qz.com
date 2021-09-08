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
Mg 10.51884 0.93935 0.79758
C 10.62390 0.46717 -2.75841
C 9.91589 -2.40650 1.07386
C 9.97404 1.39376 3.89017
C 9.99043 4.33049 0.19203
N 10.31890 -0.65827 -0.71896
C 10.43939 -0.73013 -2.03602
C 10.43785 -2.17470 -2.60417
C 9.64762 -2.90185 -1.41017
C 10.00680 -1.95486 -0.26698
C 8.13860 -2.88425 -1.58874
C 11.88690 -2.75329 -2.77006
C 11.97672 -4.03323 -3.60398
H 13.14800 -4.22958 -4.45370
N 10.50469 -0.35403 2.30274
C 10.26401 -1.68015 2.23385
C 10.34940 -2.17688 3.61593
C 10.49710 -1.00866 4.52596
C 10.30593 0.08073 3.58959
C 10.36891 -3.63901 3.95535
C 10.75359 -0.89155 5.99793
O 10.60721 0.21913 6.62066
C 11.09052 -2.13536 6.70568
N 9.84050 2.57697 1.88405
C 9.74361 2.52541 3.22415
C 9.50096 3.89733 3.80669
C 9.69682 4.92021 2.59304
C 9.94185 3.91175 1.44920
C 8.12772 4.11742 4.44084
C 10.72274 6.06818 2.79875
C 12.12192 5.50154 2.87994
N 10.51898 2.22819 -0.93648
C 10.24808 3.61368 -0.96952
C 10.25421 4.02479 -2.38846
C 10.50747 2.83611 -3.08945
C 10.62914 1.75653 -2.17336
C 9.99669 5.33730 -2.87540
C 10.68901 2.24726 -4.39673
O 10.77875 2.68651 -5.54248
C 10.64688 0.66527 -4.24867
C 11.81044 0.09604 -5.01327
O 12.97781 0.32599 -4.65137
O 11.35391 -0.79770 -5.90948
C 12.47196 -1.39800 -6.68498
H 9.49892 -3.39650 1.26072
H 9.88440 1.62575 4.95395
H 9.93108 5.38810 -0.07370
H 9.89780 -2.29595 -3.53552
H 10.04502 -3.89569 -1.22117
H 7.83979 -2.28029 -2.43259
H 7.73843 -2.47534 -0.67021
H 7.75501 -3.90634 -1.65591
H 12.23286 -3.13567 -1.80113
H 12.54925 -2.01210 -3.17324
H 11.05594 -4.22251 -4.16904
H 12.05412 -4.88543 -2.92444
H 10.07880 -4.25546 3.09556
H 9.57453 -3.80723 4.69055
H 11.32089 -3.95921 4.38415
H 10.28393 -2.85338 6.57796
H 11.32592 -1.83499 7.72854
H 11.97927 -2.57401 6.23844
H 10.22692 4.12393 4.57789
H 8.76318 5.44901 2.35114
H 8.18164 4.43838 5.49827
H 7.55904 3.20151 4.34999
H 7.49164 4.83621 3.92709
H 10.53780 6.57565 3.74871
H 10.68958 6.82194 2.01918
H 12.20264 5.19240 3.92637
H 12.85879 6.26403 2.63567
H 12.21844 4.60133 2.26190
H 10.73320 5.53471 -3.64291
H 10.04623 6.07450 -2.07674
H 8.97874 5.20900 -3.25302
H 9.68995 0.32444 -4.64579
H 12.25170 -1.23168 -7.74400
H 12.54213 -2.44802 -6.34395
H 13.44809 -0.93556 -6.59683

