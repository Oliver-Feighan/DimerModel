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
Mg 8.93053 0.80864 0.68762
C 6.69493 -1.76572 1.80737
C 6.23580 2.89739 0.58666
C 10.78045 2.96545 -0.72425
C 11.16400 -1.71279 -0.20664
N 6.79329 0.49048 1.16014
C 6.07073 -0.51260 1.63561
C 4.61592 -0.13162 2.02023
C 4.44469 1.15243 1.07114
C 5.90987 1.56766 0.95500
C 3.93374 0.80986 -0.31844
C 4.47920 0.25182 3.53558
C 3.04061 0.33482 4.05067
H 2.76648 -0.08857 5.42110
N 8.69400 2.76872 0.48584
C 7.53071 3.45282 0.49150
C 7.89282 4.86394 0.28712
C 9.35040 4.93673 -0.00422
C 9.66776 3.52709 -0.11514
C 6.91206 5.98808 0.45395
C 10.29880 6.08533 -0.16989
O 11.47404 5.92355 -0.65515
C 9.77835 7.41884 0.16533
N 10.62088 0.64772 -0.51202
C 11.26394 1.74392 -0.95106
C 12.59883 1.38616 -1.55971
C 12.82615 -0.15963 -1.21846
C 11.48381 -0.46419 -0.51793
C 12.71983 1.59706 -3.06877
C 14.14380 -0.53208 -0.48514
C 14.13364 0.03703 0.91536
N 9.03784 -1.33474 0.93753
C 10.04855 -2.19928 0.46305
C 9.62463 -3.57893 0.77814
C 8.37605 -3.42150 1.39843
C 8.03594 -2.04256 1.44708
C 10.33981 -4.76470 0.44863
C 7.26358 -4.10904 2.01357
O 7.03107 -5.27474 2.33109
C 6.06798 -3.07703 2.19210
C 5.53934 -3.21649 3.59339
O 6.24780 -2.91279 4.56933
O 4.22008 -3.47606 3.54285
C 3.64788 -3.61195 4.90884
H 5.42629 3.57019 0.30219
H 11.50504 3.67563 -1.12933
H 11.84979 -2.54783 -0.36494
H 3.87054 -0.87823 1.77373
H 3.86538 1.93098 1.56065
H 3.88224 -0.25706 -0.47708
H 4.63383 1.26167 -1.00894
H 2.97321 1.30104 -0.49859
H 4.78738 1.29670 3.66974
H 5.04523 -0.41997 4.15115
H 2.32699 -0.12515 3.35631
H 2.74366 1.38591 4.08292
H 5.88083 5.61770 0.50671
H 6.94191 6.57802 -0.46852
H 7.15569 6.63799 1.29711
H 8.90770 7.63194 -0.45030
H 10.62841 8.09705 0.06886
H 9.43351 7.40853 1.20533
H 13.39619 1.96599 -1.11125
H 12.85346 -0.77790 -2.12776
H 13.55295 2.26845 -3.35023
H 11.77869 1.98442 -3.43583
H 12.83915 0.68697 -3.65416
H 15.00057 -0.08527 -0.99550
H 14.32168 -1.60101 -0.43278
H 14.46581 1.06859 0.76497
H 14.82231 -0.50994 1.55613
H 13.11461 0.07420 1.31769
H 10.29959 -5.40481 1.31991
H 11.36635 -4.54366 0.16378
H 9.74671 -5.12136 -0.39771
H 5.31233 -3.31842 1.44356
H 3.15530 -4.58821 4.95323
H 3.01124 -2.71937 5.05622
H 4.34193 -3.68737 5.73774

