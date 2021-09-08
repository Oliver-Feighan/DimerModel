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
Mg 15.77258 1.42639 1.20254
C 15.39723 0.40455 4.62219
C 14.54473 4.46953 2.13328
C 15.66624 2.21222 -1.86969
C 15.87234 -1.95118 0.34944
N 15.09438 2.21410 3.15493
C 15.05545 1.75423 4.39645
C 14.70327 2.82170 5.46698
C 13.91666 3.85703 4.52453
C 14.57933 3.52420 3.18933
C 12.43151 3.55665 4.40866
C 15.98051 3.47741 6.10002
C 15.72166 4.30011 7.36411
H 16.73135 4.29152 8.41906
N 15.66175 3.20826 0.33581
C 15.15422 4.34083 0.86588
C 15.28670 5.36839 -0.17850
C 15.76266 4.71085 -1.42582
C 15.69394 3.31785 -1.03221
C 15.04927 6.82773 0.08105
C 16.19444 5.24637 -2.75735
O 16.34154 4.48140 -3.77520
C 16.34952 6.70487 -2.85838
N 15.55721 0.30300 -0.53353
C 15.59688 0.88680 -1.74412
C 15.69827 -0.14347 -2.84350
C 15.96323 -1.53640 -2.10383
C 15.87654 -1.06396 -0.63595
C 14.47345 -0.27313 -3.74870
C 17.21445 -2.34240 -2.54824
C 18.47214 -1.60286 -2.15240
N 15.84420 -0.44744 2.27572
C 15.85448 -1.74696 1.72325
C 15.78998 -2.69975 2.85040
C 15.72290 -1.88174 3.98832
C 15.72372 -0.51565 3.59682
C 15.74978 -4.11825 2.73933
C 15.64086 -1.86334 5.43110
O 15.69286 -2.71749 6.31511
C 15.29877 -0.38332 5.89898
C 16.23540 -0.02212 7.01915
O 17.45722 0.08282 6.81241
O 15.51521 0.35065 8.09281
C 16.40012 0.72900 9.22666
H 13.96044 5.37980 2.27072
H 15.74003 2.42803 -2.93808
H 15.99786 -3.02141 0.17188
H 14.05276 2.47249 6.25987
H 14.12485 4.88556 4.80763
H 12.17066 2.62390 4.88623
H 12.22271 3.50994 3.34788
H 11.84592 4.39607 4.79441
H 16.34515 4.26858 5.43232
H 16.73059 2.73666 6.29827
H 14.72623 4.10950 7.78332
H 15.69946 5.35795 7.09118
H 14.55063 6.98898 1.04477
H 14.32120 7.17038 -0.66217
H 15.95831 7.42563 -0.01295
H 15.40644 7.18632 -2.61086
H 16.74912 6.88773 -3.85776
H 17.07674 7.03225 -2.10699
H 16.53366 0.07101 -3.49865
H 15.13398 -2.24280 -2.25609
H 14.70417 -0.12000 -4.81986
H 13.72636 0.43501 -3.41563
H 13.94274 -1.22080 -3.67468
H 17.23837 -2.43487 -3.63683
H 17.24923 -3.34650 -2.13890
H 18.60151 -0.88291 -2.96606
H 19.31616 -2.28704 -2.09284
H 18.31916 -1.03253 -1.22883
H 16.42387 -4.50857 3.49017
H 16.03178 -4.44705 1.74126
H 14.69215 -4.30033 2.94795
H 14.25587 -0.37343 6.21809
H 16.10488 0.11299 10.08169
H 16.29477 1.82485 9.33507
H 17.45281 0.48523 9.14383

