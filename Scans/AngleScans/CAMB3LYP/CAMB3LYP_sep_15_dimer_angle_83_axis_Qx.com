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
Mg 7.87916 0.71765 0.60460
C 10.12127 0.71024 3.40674
C 5.35179 0.11985 2.81586
C 5.94956 0.27093 -1.87421
C 10.65069 0.17434 -1.43724
N 7.86402 0.46827 2.80229
C 8.76752 0.54718 3.76760
C 8.17815 0.52611 5.20348
C 6.79601 -0.22701 4.88507
C 6.63397 0.16905 3.41894
C 6.90399 -1.74087 4.96175
C 7.92936 1.96855 5.76870
C 7.64770 2.03062 7.27156
H 8.15080 3.17436 8.02734
N 5.89772 0.74938 0.49743
C 5.02472 0.50512 1.49717
C 3.68546 0.63338 0.90212
C 3.84010 0.81098 -0.56732
C 5.26538 0.59170 -0.71079
C 2.42937 0.66188 1.72362
C 2.86438 1.11379 -1.66389
O 3.18013 0.98834 -2.89977
C 1.50046 1.47138 -1.24755
N 8.21309 0.07570 -1.34443
C 7.20810 0.02096 -2.23602
C 7.73147 -0.20166 -3.63492
C 9.31935 -0.04291 -3.52989
C 9.45260 0.16486 -2.00526
C 7.39594 -1.55402 -4.26335
C 9.98757 0.99030 -4.47786
C 9.56763 2.39069 -4.09323
N 10.02233 0.66471 0.87302
C 10.99876 0.39450 -0.11068
C 12.30740 0.35646 0.57375
C 11.99652 0.58618 1.92248
C 10.59114 0.73536 2.07132
C 13.56122 0.08525 -0.04308
C 12.53662 0.72582 3.25570
O 13.66940 0.77997 3.73279
C 11.33437 0.68655 4.29467
C 11.51791 1.82483 5.26067
O 11.44140 3.00327 4.87076
O 11.53641 1.33931 6.51538
C 11.70533 2.43172 7.51033
H 4.52288 -0.29217 3.39208
H 5.34033 0.21507 -2.77931
H 11.57285 0.10788 -2.01851
H 8.75461 -0.04312 5.92282
H 5.98075 0.17500 5.48098
H 7.92431 -2.06540 5.10238
H 6.51518 -2.11148 4.02232
H 6.23695 -2.12756 5.73754
H 6.97310 2.34439 5.38265
H 8.74764 2.61874 5.52726
H 7.90439 1.09230 7.77813
H 6.56937 2.12839 7.41877
H 2.61707 0.34416 2.75669
H 1.76496 -0.10917 1.31878
H 1.91924 1.62651 1.67983
H 1.07821 0.66035 -0.65900
H 0.97520 1.73920 -2.16638
H 1.55370 2.34450 -0.58765
H 7.34756 0.54846 -4.31539
H 9.83930 -0.98309 -3.76555
H 6.85600 -1.46648 -5.22505
H 6.81540 -2.12632 -3.55203
H 8.25007 -2.20512 -4.44150
H 9.64898 0.83591 -5.50532
H 11.07096 0.93349 -4.48383
H 8.60029 2.50393 -4.59175
H 10.28759 3.12038 -4.45825
H 9.39301 2.46668 -3.01370
H 14.26788 0.79749 0.36171
H 13.49590 0.16075 -1.12651
H 13.72268 -0.94421 0.28767
H 11.36422 -0.28278 4.79374
H 12.57989 2.17843 8.11751
H 10.73431 2.51137 8.03428
H 11.98303 3.41066 7.13751

