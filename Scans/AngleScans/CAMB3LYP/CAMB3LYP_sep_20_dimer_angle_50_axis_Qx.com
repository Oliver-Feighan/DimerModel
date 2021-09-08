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
Mg 10.51137 0.95493 0.80453
C 12.21294 -0.34116 3.68617
C 7.91807 1.74270 2.87547
C 8.80900 1.69057 -1.77000
C 12.68888 -0.93084 -1.15711
N 10.22157 0.68644 2.98080
C 10.96177 0.23813 3.98354
C 10.36017 0.49373 5.39148
C 8.81386 0.61157 4.97448
C 8.98505 1.07677 3.52986
C 8.09002 -0.72457 4.96341
C 10.88524 1.82594 6.03298
C 10.58294 1.98419 7.52474
H 11.56915 2.65479 8.36745
N 8.86624 2.04990 0.62246
C 7.93484 2.28305 1.57084
C 6.91418 3.12897 0.93300
C 7.23551 3.23997 -0.51584
C 8.32820 2.29344 -0.61678
C 5.81764 3.80328 1.70522
C 6.64723 4.05275 -1.62918
O 6.92689 3.81470 -2.85717
C 5.66237 5.07461 -1.24531
N 10.57668 0.29340 -1.16536
C 9.75910 0.81440 -2.09706
C 10.17233 0.38780 -3.48545
C 11.58797 -0.33506 -3.30929
C 11.71169 -0.27770 -1.77085
C 9.20737 -0.55253 -4.20727
C 12.76550 0.20516 -4.16618
C 13.13570 1.59940 -3.71401
N 12.27087 -0.24973 1.15242
C 13.01314 -0.97251 0.19289
C 14.05042 -1.72864 0.92415
C 13.82326 -1.40875 2.27121
C 12.70948 -0.53215 2.37402
C 15.00186 -2.61245 0.34127
C 14.26578 -1.62169 3.63045
O 15.21779 -2.19933 4.15365
C 13.16412 -1.04001 4.61756
C 13.86447 -0.20861 5.65718
O 14.45587 0.83741 5.33673
O 13.53833 -0.66566 6.87992
C 14.19993 0.13406 7.94507
H 6.96176 1.82358 3.39298
H 8.32508 1.99832 -2.69984
H 13.46813 -1.46493 -1.70502
H 10.49426 -0.31757 6.09688
H 8.30320 1.37052 5.56149
H 8.76672 -1.55072 5.12371
H 7.62566 -0.79950 3.98883
H 7.27056 -0.71538 5.68779
H 10.30598 2.66833 5.63350
H 11.93812 1.94137 5.86341
H 10.26414 1.04015 7.98312
H 9.71726 2.64172 7.63554
H 5.73828 3.40340 2.72365
H 4.87197 3.52284 1.22885
H 5.90677 4.89172 1.70016
H 4.83437 4.60028 -0.72387
H 5.42322 5.61041 -2.16596
H 6.13112 5.76183 -0.53214
H 10.29466 1.24687 -4.13363
H 11.53845 -1.39972 -3.58109
H 8.86222 -0.15946 -5.18218
H 8.36570 -0.74430 -3.55507
H 9.59022 -1.55493 -4.39118
H 12.46486 0.28812 -5.21350
H 13.64812 -0.42484 -4.13342
H 12.41400 2.22983 -4.24201
H 14.15636 1.83839 -4.00585
H 12.95870 1.72464 -2.63947
H 15.95165 -2.40429 0.81567
H 15.05803 -2.48096 -0.73731
H 14.56557 -3.57680 0.61505
H 12.63813 -1.88803 5.05775
H 14.76147 -0.56782 8.56935
H 13.39033 0.70725 8.43471
H 14.98188 0.82106 7.64336

