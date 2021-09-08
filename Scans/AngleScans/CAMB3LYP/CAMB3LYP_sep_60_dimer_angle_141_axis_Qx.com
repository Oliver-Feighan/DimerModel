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
Mg 31.52396 2.84036 2.40133
C 32.83903 5.04412 4.90993
C 30.78319 0.66017 4.91780
C 30.77291 0.69505 0.18751
C 33.36751 4.63850 0.05320
N 31.82056 2.95136 4.59038
C 32.27418 3.86384 5.43653
C 32.04041 3.52565 6.93337
C 31.93112 1.92865 6.80485
C 31.44747 1.83229 5.35935
C 33.27328 1.22337 6.90746
C 30.71132 4.14968 7.48665
C 30.57328 4.12164 9.01043
H 29.90409 5.23987 9.66938
N 30.44127 1.18247 2.53458
C 30.22740 0.43683 3.63892
C 29.38294 -0.68791 3.20761
C 29.25211 -0.63511 1.72612
C 30.18785 0.42728 1.41656
C 28.72733 -1.63071 4.17442
C 28.43143 -1.42055 0.74850
O 28.65257 -1.36648 -0.51285
C 27.42276 -2.32611 1.31776
N 32.16167 2.55194 0.44365
C 31.63687 1.57577 -0.31766
C 32.04350 1.73336 -1.76342
C 32.75604 3.16210 -1.85602
C 32.71571 3.56212 -0.36480
C 32.98362 0.65965 -2.31088
C 32.19556 4.16141 -2.90479
C 30.80357 4.59844 -2.50884
N 32.71728 4.64146 2.40900
C 33.42225 5.20224 1.32151
C 34.17789 6.35966 1.84285
C 33.87585 6.37826 3.21286
C 33.01025 5.29594 3.52717
C 35.04647 7.19541 1.08575
C 34.10092 7.06108 4.46658
O 34.67638 8.09581 4.80111
C 33.54049 6.15277 5.64442
C 32.71533 7.02452 6.55088
O 31.66051 7.54135 6.14242
O 33.18956 6.92824 7.80635
C 32.39670 7.76672 8.74451
H 30.71673 -0.18706 5.60093
H 30.45845 0.04871 -0.63519
H 33.88832 5.30901 -0.63372
H 32.85876 3.79052 7.59210
H 31.18360 1.52792 7.48455
H 34.09534 1.92321 6.93170
H 33.34077 0.59062 6.03218
H 33.27974 0.54867 7.76831
H 29.86936 3.50224 7.20977
H 30.58476 5.15365 7.13077
H 31.52537 3.89722 9.50651
H 29.92465 3.28616 9.28479
H 29.13984 -1.52168 5.18504
H 29.01035 -2.64521 3.87351
H 27.63817 -1.55095 4.16762
H 27.91038 -3.04288 1.97422
H 26.87826 -2.73154 0.46276
H 26.73993 -1.74033 1.94323
H 31.17582 1.73072 -2.41168
H 33.81782 3.07100 -2.12834
H 32.58212 0.14107 -3.20191
H 33.19036 -0.04868 -1.51963
H 33.98041 1.00928 -2.57424
H 32.10290 3.67548 -3.87918
H 32.81820 5.03938 -3.04068
H 30.17326 3.78900 -2.88910
H 30.55229 5.54779 -2.97750
H 30.69252 4.61825 -1.41846
H 34.83566 8.21408 1.38303
H 34.90181 7.05441 0.01665
H 36.01773 6.82215 1.42130
H 34.39818 5.72069 6.16139
H 33.10097 8.43656 9.24751
H 31.83637 7.05556 9.38014
H 31.69942 8.47665 8.31525

