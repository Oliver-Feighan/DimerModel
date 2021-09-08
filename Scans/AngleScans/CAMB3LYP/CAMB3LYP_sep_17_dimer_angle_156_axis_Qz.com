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
Mg 8.94232 0.79850 0.67509
C 8.92026 1.67432 -2.80507
C 8.42316 -2.43219 -0.28811
C 8.50756 0.06186 3.73055
C 8.31754 4.15371 1.38074
N 8.71868 -0.13344 -1.31833
C 8.78978 0.28945 -2.57150
C 8.79762 -0.84308 -3.63304
C 8.07001 -1.98589 -2.77081
C 8.45244 -1.51524 -1.36906
C 6.55520 -1.95561 -2.88812
C 10.25139 -1.26955 -4.04134
C 10.33668 -2.14805 -5.29143
H 11.47827 -1.97693 -6.18582
N 9.01432 -0.95831 1.59542
C 8.79996 -2.17308 1.04802
C 8.94938 -3.14105 2.14569
C 9.10674 -2.38622 3.41856
C 8.85600 -1.03565 2.95684
C 9.01367 -4.62373 1.92000
C 9.41725 -2.81122 4.82182
O 9.27097 -2.01405 5.81478
C 9.80813 -4.21595 5.01007
N 8.27107 1.89591 2.30833
C 8.22712 1.35068 3.53664
C 7.97746 2.40203 4.59150
C 8.10411 3.80634 3.83688
C 8.32661 3.29963 2.39496
C 6.62527 2.32547 5.30009
C 9.11208 4.83220 4.42364
C 10.52530 4.32402 4.25063
N 8.84759 2.63497 -0.45912
C 8.54564 3.92488 0.02982
C 8.48809 4.83004 -1.13629
C 8.73982 3.99294 -2.23373
C 8.92012 2.65657 -1.78529
C 8.18357 6.21993 -1.09627
C 8.88348 3.93408 -3.67060
O 8.91939 4.76756 -4.57490
C 8.88140 2.40847 -4.11660
C 10.02665 2.20157 -5.06955
O 11.20186 2.32194 -4.68096
O 9.55536 1.68601 -6.21953
C 10.65540 1.45269 -7.19282
H 8.03527 -3.43506 -0.46884
H 8.45405 -0.11786 4.80680
H 8.22506 5.23208 1.52648
H 8.22478 -0.63097 -4.52782
H 8.49587 -2.96510 -2.97367
H 6.21101 -1.09378 -3.44032
H 6.18201 -1.92812 -1.87273
H 6.19155 -2.89352 -3.31757
H 10.64270 -1.97004 -3.29238
H 10.88146 -0.40963 -4.16031
H 9.39909 -2.14726 -5.86053
H 10.45871 -3.18753 -4.97740
H 8.70403 -4.88941 0.90170
H 8.25209 -5.07834 2.56280
H 9.98821 -5.04650 2.17339
H 9.01297 -4.86362 4.64855
H 10.07627 -4.30588 6.06453
H 10.68749 -4.42062 4.38922
H 8.72756 2.35321 5.37129
H 7.15058 4.35462 3.83367
H 6.71300 2.23568 6.39932
H 6.07350 1.48857 4.89304
H 5.95440 3.16069 5.10631
H 8.95298 4.94707 5.49865
H 9.03251 5.81863 3.97909
H 10.65303 3.65389 5.10602
H 11.23549 5.14773 4.28504
H 10.61739 3.71893 3.34121
H 8.88547 6.71157 -1.75668
H 8.24792 6.61204 -0.08342
H 7.15483 6.20492 -1.46605
H 7.91748 2.20534 -4.58475
H 10.39086 1.99001 -8.10879
H 10.76142 0.35400 -7.26611
H 11.62394 1.88336 -6.96722

