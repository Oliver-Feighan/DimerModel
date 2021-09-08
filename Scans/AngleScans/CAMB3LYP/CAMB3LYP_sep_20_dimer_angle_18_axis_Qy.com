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
Mg 10.51917 0.95107 0.80294
C 11.11153 -1.20554 3.60957
C 10.02506 3.47851 3.03957
C 9.55309 2.76676 -1.61313
C 9.93427 -1.91177 -1.09648
N 10.58760 1.00550 3.01305
C 10.88350 0.13890 3.97002
C 11.01707 0.76368 5.38472
C 10.10738 2.06760 5.15860
C 10.27298 2.22224 3.64819
C 8.63564 1.83364 5.45626
C 12.49818 1.14950 5.73034
C 12.75039 1.47671 7.20380
H 14.03060 1.09353 7.79263
N 10.36071 2.92433 0.66266
C 10.18669 3.80065 1.67412
C 10.10983 5.12907 1.04677
C 10.08749 4.94964 -0.43039
C 9.97614 3.50748 -0.51910
C 10.14537 6.40519 1.83651
C 10.14265 5.91698 -1.57382
O 9.86930 5.55780 -2.77344
C 10.43491 7.31728 -1.23437
N 9.63589 0.51029 -1.02678
C 9.36547 1.48205 -1.91586
C 8.99120 0.90334 -3.25952
C 9.30248 -0.66183 -3.15495
C 9.73794 -0.73491 -1.67494
C 7.53839 1.10541 -3.68951
C 10.24359 -1.26578 -4.23301
C 11.64226 -0.71803 -4.06249
N 10.69411 -1.18066 1.10823
C 10.37240 -2.20322 0.18901
C 10.54834 -3.49101 0.89120
C 10.94110 -3.12472 2.18737
C 10.98587 -1.70819 2.29189
C 10.30275 -4.77985 0.33920
C 11.31772 -3.61409 3.49403
O 11.53928 -4.72461 3.97507
C 11.32243 -2.39136 4.50956
C 12.59931 -2.45092 5.30236
O 13.69820 -2.28541 4.74409
O 12.30253 -2.47601 6.61448
C 13.53708 -2.52410 7.44237
H 9.62874 4.28500 3.65712
H 9.31445 3.34389 -2.50949
H 9.86643 -2.85191 -1.64803
H 10.60999 0.16101 6.18774
H 10.51748 2.93015 5.67743
H 8.42573 0.79372 5.65785
H 8.10122 2.16356 4.57501
H 8.30719 2.48656 6.26998
H 12.72908 2.12323 5.27946
H 13.17642 0.38304 5.40917
H 11.92081 1.15534 7.84525
H 12.77303 2.56308 7.31908
H 9.99705 6.21972 2.90749
H 9.26873 6.99034 1.53822
H 11.04504 6.99398 1.64525
H 9.68283 7.68520 -0.54050
H 10.52074 7.83787 -2.19014
H 11.39500 7.35655 -0.70777
H 9.59910 1.33066 -4.04752
H 8.38764 -1.26602 -3.24355
H 7.43906 1.62415 -4.66169
H 7.02532 1.65465 -2.91128
H 6.94616 0.19451 -3.75826
H 9.91408 -0.96995 -5.23204
H 10.28206 -2.34983 -4.21445
H 11.59718 0.24132 -4.58660
H 12.37246 -1.38323 -4.51895
H 11.85683 -0.50793 -3.00823
H 11.12512 -5.41207 0.64663
H 10.21584 -4.73823 -0.74454
H 9.34969 -5.01895 0.81855
H 10.44123 -2.48963 5.14471
H 13.45155 -3.39995 8.09295
H 13.60539 -1.53576 7.93457
H 14.47294 -2.72786 6.93532

