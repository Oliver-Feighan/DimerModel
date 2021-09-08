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
Mg 7.88535 0.71240 0.58882
C 10.43436 -1.68416 -0.21002
C 10.18483 3.14888 -0.05193
C 5.80867 2.83016 1.71571
C 6.18878 -1.84813 2.23535
N 10.00308 0.61910 -0.04268
C 10.84737 -0.34569 -0.37537
C 12.18235 0.15203 -0.99154
C 12.20452 1.62733 -0.35760
C 10.70326 1.83439 -0.16779
C 12.87105 1.68610 1.00682
C 12.13553 0.20110 -2.55915
C 13.49534 0.38014 -3.23785
H 13.73196 -0.29460 -4.51112
N 7.79944 2.68044 0.34797
C 8.82538 3.51166 0.06908
C 8.24501 4.86227 0.01267
C 6.82292 4.78009 0.44352
C 6.75962 3.40028 0.88193
C 9.00439 6.05371 -0.49449
C 5.70829 5.78164 0.47006
O 4.61835 5.55937 1.10681
C 5.96693 7.06808 -0.19302
N 6.34381 0.57457 1.97699
C 5.55899 1.63271 2.24589
C 4.35500 1.22560 3.06141
C 4.36876 -0.37348 3.07900
C 5.68581 -0.62367 2.31209
C 4.31559 1.74977 4.49671
C 3.08161 -1.09376 2.59205
C 2.88763 -0.85902 1.11129
N 8.12574 -1.42366 0.80685
C 7.31626 -2.30046 1.56187
C 7.94354 -3.63641 1.49723
C 9.09595 -3.43772 0.72201
C 9.19160 -2.07008 0.34798
C 7.46799 -4.81096 2.14562
C 10.25890 -4.07183 0.14395
O 10.66152 -5.23091 0.05357
C 11.24414 -2.93963 -0.37932
C 11.68010 -3.30989 -1.77050
O 10.85777 -3.34263 -2.70286
O 13.02424 -3.35140 -1.81465
C 13.50508 -3.70379 -3.17716
H 10.88731 3.98158 -0.00705
H 5.00687 3.49680 2.04158
H 5.66931 -2.71906 2.64070
H 13.06060 -0.40013 -0.67916
H 12.60367 2.35351 -1.06095
H 13.11542 0.70212 1.37855
H 12.15913 2.17176 1.66131
H 13.74485 2.34303 0.97151
H 11.64462 1.13102 -2.87378
H 11.64676 -0.66868 -2.95314
H 14.32781 0.19990 -2.54688
H 13.60564 1.42935 -3.52248
H 10.07653 5.83926 -0.58357
H 8.94515 6.82463 0.28142
H 8.58992 6.45360 -1.42238
H 6.83307 7.54170 0.26299
H 5.02371 7.61582 -0.14619
H 6.22796 6.87592 -1.23971
H 3.43864 1.56402 2.59357
H 4.51671 -0.76853 4.09482
H 3.40431 2.33506 4.72276
H 5.20231 2.34584 4.66687
H 4.39807 0.98718 5.26927
H 2.20284 -0.67682 3.09013
H 3.08484 -2.16137 2.78491
H 2.39716 0.11846 1.07856
H 2.25518 -1.63258 0.68042
H 3.85166 -0.76143 0.59868
H 7.54988 -5.61714 1.42877
H 6.44364 -4.68745 2.49089
H 8.17611 -4.87558 2.97608
H 12.08484 -2.89096 0.31381
H 14.15183 -4.58022 -3.07098
H 13.96731 -2.78234 -3.57844
H 12.77270 -4.06638 -3.88892

