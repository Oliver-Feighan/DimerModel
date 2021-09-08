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
Mg 13.14629 1.18410 0.99091
C 12.69532 4.63306 0.10756
C 12.58034 0.42562 -2.28613
C 13.09076 -1.88006 1.81268
C 12.54365 2.09678 4.29992
N 12.70665 2.45263 -0.76687
C 12.62061 3.75477 -0.99355
C 12.52501 4.14151 -2.49384
C 11.92715 2.76277 -3.05995
C 12.46346 1.80827 -1.99507
C 10.40880 2.70775 -3.02961
C 13.92833 4.45906 -3.11984
C 13.88231 5.14002 -4.48949
H 14.90645 6.13074 -4.80899
N 13.36120 -0.44722 -0.11868
C 13.10745 -0.57310 -1.43816
C 13.40472 -1.97438 -1.77314
C 13.69676 -2.70661 -0.51089
C 13.36666 -1.68896 0.46661
C 13.47095 -2.47681 -3.18615
C 14.17945 -4.09735 -0.22970
O 14.13608 -4.59979 0.94864
C 14.61737 -4.88650 -1.39029
N 12.65164 0.20698 2.75809
C 12.76408 -1.12868 2.86436
C 12.62042 -1.57570 4.29961
C 12.62926 -0.23429 5.17032
C 12.68946 0.80661 4.03089
C 11.36336 -2.38324 4.62219
C 13.67925 -0.14278 6.31129
C 15.07177 -0.07388 5.72690
N 12.88195 3.03727 2.06976
C 12.61471 3.19139 3.44777
C 12.40179 4.63283 3.69141
C 12.53826 5.21930 2.42416
C 12.79658 4.21200 1.45559
C 12.07692 5.22569 4.94404
C 12.51199 6.46026 1.68387
O 12.42424 7.64988 1.98527
C 12.48710 6.12193 0.13104
C 13.51536 6.98398 -0.54871
O 14.72571 6.83577 -0.30438
O 12.92150 7.68749 -1.52992
C 13.90298 8.54958 -2.24077
H 12.19355 0.06135 -3.23829
H 13.16853 -2.91008 2.16834
H 12.44791 2.47307 5.32066
H 11.84624 4.95775 -2.71017
H 12.34515 2.52131 -4.03378
H 9.98482 3.56081 -2.52102
H 10.15793 1.79328 -2.50820
H 10.01536 2.59404 -4.04381
H 14.41921 3.51668 -3.39494
H 14.52289 5.04396 -2.44528
H 12.88413 5.53028 -4.72245
H 14.06104 4.38477 -5.25869
H 13.04824 -1.75071 -3.89148
H 12.79987 -3.34018 -3.25031
H 14.47674 -2.78870 -3.47558
H 13.79788 -4.97005 -2.10016
H 15.01014 -5.81828 -0.97870
H 15.42092 -4.34306 -1.89979
H 13.45834 -2.19378 4.59786
H 11.67147 -0.07833 5.68815
H 11.58231 -3.37443 5.06227
H 10.78383 -2.48967 3.71480
H 10.65800 -1.89237 5.29056
H 13.64627 -1.04332 6.92943
H 13.52828 0.70580 6.97017
H 15.31297 -1.12448 5.53942
H 15.76472 0.36911 6.43937
H 15.06745 0.44111 4.75921
H 12.68592 6.11529 5.03474
H 12.25296 4.53677 5.76755
H 11.01214 5.42432 4.79551
H 11.47872 6.32401 -0.23216
H 13.52140 9.57420 -2.19345
H 14.02114 8.10198 -3.24548
H 14.88278 8.66789 -1.79305

