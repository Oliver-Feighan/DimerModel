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
Mg 10.51895 0.94072 0.79446
C 10.47260 2.03384 -2.62344
C 10.00908 -2.22543 -0.36749
C 10.10539 0.01097 3.79978
C 9.88285 4.24168 1.71353
N 10.28801 0.13506 -1.25245
C 10.34986 0.63639 -2.47687
C 10.35674 -0.42699 -3.60768
C 9.63944 -1.62543 -2.81516
C 10.02784 -1.24211 -1.38867
C 8.12384 -1.59535 -2.92196
C 11.81005 -0.81969 -4.05005
C 11.89211 -1.61729 -5.35346
H 13.02768 -1.38452 -6.24158
N 10.60440 -0.87018 1.60186
C 10.39246 -2.04912 0.98020
C 10.55272 -3.08354 2.01387
C 10.71401 -2.40959 3.33088
C 10.45437 -1.03389 2.95656
C 10.62252 -4.54873 1.69486
C 11.03464 -2.92057 4.70284
O 10.89047 -2.18825 5.74484
C 11.43307 -4.33241 4.80006
N 9.85218 1.92975 2.49727
C 9.81789 1.30804 3.68901
C 9.56954 2.28962 4.80937
C 9.68532 3.73928 4.14404
C 9.90175 3.32549 2.67185
C 8.22186 2.16188 5.51918
C 10.69195 4.73114 4.78869
C 12.10647 4.24189 4.57621
N 10.40916 2.84448 -0.22126
C 10.10413 4.09952 0.34963
C 10.03562 5.07602 -0.75680
C 10.28481 4.31096 -1.90618
C 10.47388 2.94991 -1.54382
C 9.72494 6.45910 -0.62762
C 10.42037 4.34341 -3.34467
O 10.44718 5.23236 -4.19484
C 10.42273 2.84893 -3.88588
C 11.56335 2.70814 -4.85628
O 12.74024 2.80963 -4.46737
O 11.08775 2.26368 -6.03383
C 12.18317 2.09759 -7.02595
H 9.62477 -3.21685 -0.60890
H 10.05897 -0.23643 4.86285
H 9.78625 5.30826 1.92742
H 9.77773 -0.16180 -4.48415
H 10.06861 -2.58781 -3.08165
H 7.77248 -0.70218 -3.41686
H 7.75645 -1.63371 -1.90481
H 7.76202 -2.50614 -3.40764
H 12.20894 -1.56400 -3.34888
H 12.43546 0.04914 -4.11809
H 10.95124 -1.58531 -5.91620
H 12.02076 -2.67386 -5.10622
H 10.30819 -4.75129 0.66357
H 9.86680 -5.04669 2.31194
H 11.60046 -4.98177 1.91574
H 10.63881 -4.95997 4.40284
H 11.70776 -4.48724 5.84526
H 12.30973 -4.49321 4.16270
H 10.32438 2.19552 5.58039
H 8.72927 4.28193 4.18064
H 8.31640 2.00348 6.61008
H 7.67159 1.34955 5.06327
H 7.54603 3.00430 5.38212
H 10.53859 4.77729 5.86966
H 10.60526 5.74320 4.40762
H 12.24226 3.51986 5.38697
H 12.81304 5.06531 4.65853
H 12.19605 3.69573 3.62998
H 10.42071 6.99483 -1.25960
H 9.79338 6.78696 0.40756
H 8.69415 6.46231 -0.99192
H 9.45704 2.67090 -4.36057
H 11.91083 2.69020 -7.90477
H 12.29382 1.00623 -7.16890
H 13.15102 2.51799 -6.77902

