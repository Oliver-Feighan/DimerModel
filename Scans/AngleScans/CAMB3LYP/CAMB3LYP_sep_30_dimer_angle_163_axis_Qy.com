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
Mg 15.76173 1.42162 1.18835
C 17.08368 -1.21337 -0.85818
C 17.24370 3.62513 -0.95242
C 14.99743 3.74554 3.20892
C 15.37863 -0.93257 3.72936
N 16.97895 1.12284 -0.63410
C 17.36445 0.07801 -1.35108
C 18.02130 0.43798 -2.71069
C 18.51981 1.92169 -2.35155
C 17.49947 2.26698 -1.26879
C 19.90613 1.95686 -1.73010
C 16.98642 0.44945 -3.89005
C 17.60586 0.48826 -5.28882
H 16.93907 -0.23844 -6.36574
N 15.64394 3.38239 0.90535
C 16.29635 4.11096 -0.02461
C 15.88485 5.50604 0.19583
C 15.06429 5.56267 1.43607
C 15.22424 4.20614 1.92014
C 16.20550 6.61102 -0.76834
C 14.27761 6.66012 2.08640
O 13.83627 6.55305 3.28497
C 14.11921 7.90028 1.31283
N 15.45818 1.45889 3.24476
C 15.08242 2.58966 3.86773
C 14.65816 2.31349 5.29045
C 14.59706 0.72060 5.41904
C 15.10507 0.33352 4.01281
C 15.57204 2.87696 6.37844
C 13.26064 0.10551 5.91736
C 12.17792 0.31729 4.88381
N 15.97478 -0.72083 1.36878
C 15.79035 -1.50149 2.53095
C 16.16127 -2.88938 2.18654
C 16.56088 -2.81504 0.84367
C 16.46624 -1.47181 0.38951
C 16.14959 -3.99915 3.07782
C 17.05135 -3.56559 -0.28982
O 17.24256 -4.75799 -0.52513
C 17.53200 -2.54020 -1.40514
C 16.95867 -2.98468 -2.72281
O 15.73012 -2.96855 -2.91520
O 17.96003 -3.14723 -3.60668
C 17.44039 -3.57716 -4.93218
H 17.85476 4.39265 -1.42824
H 14.62474 4.48954 3.91677
H 15.19371 -1.74257 4.43815
H 18.86636 -0.18194 -2.98519
H 18.41449 2.59042 -3.20189
H 20.28018 0.96518 -1.52322
H 19.80299 2.52135 -0.81263
H 20.58834 2.53195 -2.36271
H 16.45682 1.41083 -3.89236
H 16.31446 -0.38333 -3.81541
H 18.67694 0.25277 -5.27198
H 17.56308 1.51551 -5.65890
H 16.96041 6.29938 -1.50072
H 16.69569 7.40445 -0.19384
H 15.31530 7.02159 -1.24964
H 15.09992 8.30657 1.07702
H 13.45353 8.53113 1.90509
H 13.64091 7.65789 0.35727
H 13.67345 2.71988 5.48616
H 15.33911 0.34099 6.13674
H 15.04733 3.54708 7.08532
H 16.39228 3.39575 5.90032
H 16.08938 2.13081 6.97912
H 12.92594 0.61234 6.82580
H 13.33098 -0.95263 6.14624
H 11.83121 1.33345 5.09428
H 11.37710 -0.40785 5.01383
H 12.59545 0.31476 3.87024
H 15.71268 -4.82815 2.53710
H 15.59031 -3.77553 3.98389
H 17.22032 -4.10472 3.27189
H 18.62263 -2.54844 -1.40958
H 17.95921 -4.50473 -5.19318
H 17.58659 -2.71165 -5.60542
H 16.40472 -3.89163 -4.98651

