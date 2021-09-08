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
Mg 31.53683 2.84596 2.40187
C 32.23621 0.80666 5.27088
C 29.52856 4.72881 4.41591
C 30.54945 4.38639 -0.19034
C 32.64678 0.20010 0.42370
N 31.02923 2.68841 4.54891
C 31.40343 1.90959 5.55278
C 30.91937 2.38342 6.94946
C 29.64738 3.24575 6.48331
C 30.10389 3.61168 5.07251
C 28.37192 2.42670 6.37523
C 31.98288 3.27694 7.67923
C 31.70897 3.51855 9.16529
H 32.84595 3.60452 10.07762
N 30.63510 4.60429 2.21829
C 29.87812 5.23108 3.14317
C 29.42896 6.48285 2.51428
C 29.84660 6.46640 1.08599
C 30.35277 5.11306 0.97492
C 28.74991 7.58053 3.28071
C 29.78702 7.49393 -0.00343
O 29.98758 7.18584 -1.23131
C 29.39518 8.85293 0.39789
N 31.38871 2.29225 0.40342
C 30.97771 3.17003 -0.52852
C 31.21360 2.63702 -1.92159
C 32.09398 1.31537 -1.73255
C 32.14151 1.26149 -0.18975
C 29.95847 2.30129 -2.72678
C 33.43314 1.24372 -2.51625
C 34.40189 2.27180 -1.97778
N 32.47631 0.93161 2.74930
C 32.83289 -0.03191 1.78051
C 33.33442 -1.21583 2.50798
C 33.21233 -0.86492 3.86092
C 32.65396 0.43736 3.96940
C 33.77487 -2.43217 1.91413
C 33.41894 -1.30422 5.22218
O 33.94372 -2.28476 5.74837
C 32.67851 -0.29130 6.19796
C 33.63151 0.06784 7.30491
O 34.67055 0.70658 7.06198
O 33.05612 -0.20969 8.48908
C 33.95852 0.13973 9.61841
H 28.70125 5.24680 4.90184
H 30.32760 4.91629 -1.11957
H 33.10288 -0.62801 -0.12299
H 30.60579 1.58841 7.61535
H 29.53218 4.13972 7.09070
H 28.55713 1.37218 6.51609
H 27.98561 2.61358 5.38181
H 27.61795 2.80987 7.06880
H 31.90445 4.30546 7.30408
H 32.96915 2.87405 7.55442
H 30.94954 2.83336 9.56117
H 31.26176 4.50903 9.27916
H 28.42962 7.23966 4.27299
H 27.81510 7.80583 2.75608
H 29.35190 8.49001 3.33650
H 28.41273 8.82298 0.86304
H 29.49665 9.46384 -0.50120
H 30.09496 9.20723 1.16301
H 31.77118 3.34828 -2.51852
H 31.55397 0.41551 -2.06193
H 29.90165 2.84008 -3.69143
H 29.09252 2.52142 -2.11669
H 29.82159 1.24446 -2.94943
H 33.27019 1.49174 -3.56796
H 33.90018 0.26501 -2.48480
H 34.10426 3.18739 -2.49758
H 35.42646 1.99612 -2.21924
H 34.24566 2.43600 -0.90527
H 34.67862 -2.72264 2.43313
H 33.94919 -2.31321 0.84679
H 32.91349 -3.07305 2.12045
H 31.78509 -0.79183 6.57324
H 34.07616 -0.76375 10.22469
H 33.49804 1.01841 10.10792
H 34.99042 0.37169 9.38212

