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
Mg 10.51792 0.94285 0.80557
C 10.68747 -2.54393 -0.02669
C 9.48365 -0.00286 3.91531
C 9.92660 3.99716 1.42924
C 10.44820 1.61319 -2.61372
N 10.16951 -1.04638 1.70769
C 10.32631 -2.30160 1.31505
C 10.16056 -3.35411 2.44391
C 9.24232 -2.49390 3.44178
C 9.67876 -1.09033 3.02681
C 7.75514 -2.63380 3.16239
C 11.52890 -3.73757 3.10922
C 11.48410 -4.97942 4.00227
H 12.65005 -5.85816 4.03149
N 10.28570 1.87772 2.54069
C 9.88038 1.33663 3.70869
C 9.85281 2.44386 4.67691
C 10.11610 3.71303 3.94574
C 10.09868 3.23721 2.57703
C 9.67362 2.23063 6.15199
C 10.33326 5.12398 4.40210
O 10.30785 6.10846 3.58174
C 10.48352 5.32885 5.85018
N 10.01608 2.55745 -0.40422
C 9.86611 3.78915 0.11380
C 9.78005 4.83204 -0.97495
C 10.14756 4.06852 -2.33125
C 10.30220 2.63401 -1.78019
C 8.42456 5.52115 -1.13040
C 11.30364 4.66214 -3.18200
C 12.61535 4.51396 -2.44521
N 10.74376 -0.20683 -1.00973
C 10.65317 0.26759 -2.33668
C 10.76184 -0.90621 -3.22708
C 10.88599 -1.99331 -2.34884
C 10.83648 -1.53191 -1.00565
C 10.69124 -0.88147 -4.64853
C 11.03700 -3.42740 -2.25252
O 11.22248 -4.33749 -3.05941
C 10.78818 -3.86155 -0.74388
C 11.89590 -4.79584 -0.34069
O 13.06954 -4.39088 -0.26880
O 11.36103 -5.94863 0.10134
C 12.42006 -6.90373 0.52330
H 8.93744 -0.18153 4.84193
H 9.82976 5.07386 1.58676
H 10.53369 1.74919 -3.69393
H 9.64270 -4.25887 2.14895
H 9.50275 -2.68285 4.47997
H 7.56583 -3.19768 2.26109
H 7.37820 -1.62427 3.06327
H 7.24690 -3.06130 4.03145
H 11.78880 -2.97751 3.85722
H 12.29423 -3.85388 2.36668
H 10.56724 -5.56231 3.85205
H 11.42809 -4.65559 5.04439
H 9.33780 1.20996 6.37301
H 8.83878 2.86561 6.46777
H 10.56121 2.50166 6.72764
H 9.59694 4.96065 6.36092
H 10.71903 6.38783 5.97261
H 11.32500 4.72266 6.20382
H 10.50127 5.62318 -0.81053
H 9.29815 4.04697 -3.02973
H 8.48159 6.62222 -1.03842
H 7.74725 5.11248 -0.39231
H 7.90393 5.31102 -2.06314
H 11.15150 5.73376 -3.33245
H 11.39448 4.20884 -4.16346
H 12.61894 5.37632 -1.77184
H 13.45164 4.55248 -3.14032
H 12.61801 3.61072 -1.82408
H 11.47358 -1.53522 -5.01053
H 10.80618 0.12916 -5.03477
H 9.67923 -1.26628 -4.80027
H 9.81023 -4.34230 -0.69838
H 12.26038 -7.82764 -0.04112
H 12.34360 -6.96693 1.62508
H 13.44348 -6.66719 0.25663

