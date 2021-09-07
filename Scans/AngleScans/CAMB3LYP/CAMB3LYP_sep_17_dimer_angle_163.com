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
Mg 8.94233 0.79908 0.67430
C 8.87105 2.10503 -2.66762
C 8.44029 -2.29011 -0.68198
C 8.55055 -0.32012 3.61732
C 8.29638 4.03242 1.80301
N 8.70318 0.12279 -1.41799
C 8.75559 0.70053 -2.60871
C 8.76079 -0.28948 -3.80429
C 8.05364 -1.53900 -3.08491
C 8.44858 -1.24434 -1.63927
C 6.53733 -1.50977 -3.18123
C 10.21330 -0.64632 -4.27854
C 10.29144 -1.35984 -5.63004
H 11.42073 -1.06596 -6.50799
N 9.04081 -1.05863 1.36553
C 8.83069 -2.19712 0.67201
C 9.00173 -3.29376 1.63755
C 9.16758 -2.70330 2.99349
C 8.89943 -1.30806 2.70803
C 9.07642 -4.73559 1.22649
C 9.49854 -3.29807 4.32875
O 9.35706 -2.63355 5.41559
C 9.90403 -4.71123 4.33462
N 8.28093 1.67559 2.43970
C 8.25644 0.97991 3.59007
C 8.01010 1.88773 4.77141
C 8.11532 3.37697 4.19810
C 8.32509 3.05779 2.70159
C 6.66718 1.70886 5.47919
C 9.12111 4.33125 4.89839
C 10.53661 3.86338 4.64770
N 8.81786 2.76246 -0.21883
C 8.51038 3.97752 0.43167
C 8.43093 5.02141 -0.61064
C 8.67694 4.33150 -1.80730
C 8.87438 2.95129 -1.53246
C 8.11464 6.39206 -0.39287
C 8.80398 4.45515 -3.24162
O 8.82174 5.39599 -4.03423
C 8.81007 2.99775 -3.87591
C 9.94566 2.92401 -4.85953
O 11.12431 3.00665 -4.47151
O 9.46526 2.55223 -6.06010
C 10.55564 2.45439 -7.06675
H 8.05915 -3.26625 -0.98325
H 8.51146 -0.63420 4.66293
H 8.19613 5.08292 2.08421
H 8.17547 0.02750 -4.65910
H 8.48568 -2.48053 -3.41386
H 6.17898 -0.58896 -3.61693
H 6.17604 -1.61394 -2.16651
H 6.17688 -2.38995 -3.72132
H 10.61969 -1.43132 -3.62787
H 10.83429 0.22816 -4.29516
H 9.34717 -1.29715 -6.18445
H 10.42639 -2.42923 -5.45058
H 8.75702 -4.87436 0.18624
H 8.32660 -5.27517 1.81513
H 10.05761 -5.17682 1.41424
H 9.11037 -5.31646 3.90305
H 10.18550 -4.93020 5.36646
H 10.77770 -4.82722 3.68357
H 8.76984 1.74901 5.53078
H 7.15702 3.91149 4.27409
H 6.76879 1.48255 6.55740
H 6.11801 0.92413 4.97604
H 5.98669 2.55488 5.39920
H 8.97383 4.30848 5.98095
H 9.02754 5.36486 4.58233
H 10.68044 3.09241 5.41060
H 11.23985 4.68348 4.77783
H 10.62319 3.37835 3.66847
H 8.80424 6.96997 -0.99368
H 8.18759 6.65441 0.66051
H 7.08174 6.41308 -0.75056
H 7.84246 2.84517 -4.35553
H 10.27547 3.09981 -7.90497
H 10.67048 1.37476 -7.27878
H 11.52296 2.86321 -6.79916

