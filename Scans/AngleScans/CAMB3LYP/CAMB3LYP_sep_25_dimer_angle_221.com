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
Mg 13.14625 1.18458 0.99096
C 12.67426 4.68005 0.32916
C 12.56472 0.63118 -2.32418
C 13.10961 -1.92548 1.61835
C 12.55868 1.88407 4.35417
N 12.69054 2.55910 -0.68095
C 12.59719 3.87248 -0.82467
C 12.49107 4.35246 -2.29706
C 11.89628 3.00915 -2.94560
C 12.44318 1.99216 -1.94595
C 10.37840 2.94479 -2.91041
C 13.88924 4.71577 -2.90954
C 13.83212 5.48138 -4.23329
H 14.84980 6.49533 -4.49540
N 13.36220 -0.37254 -0.22037
C 13.10136 -0.41633 -1.54374
C 13.40312 -1.79224 -1.96797
C 13.70587 -2.60105 -0.75599
C 13.37679 -1.64863 0.28549
C 13.46344 -2.20435 -3.41017
C 14.19660 -4.00432 -0.56565
O 14.16240 -4.58017 0.57892
C 14.63138 -4.71663 -1.77605
N 12.66639 0.09568 2.69577
C 12.78560 -1.24345 2.71706
C 12.65236 -1.78067 4.12206
C 12.66009 -0.49674 5.07549
C 12.70886 0.61414 4.00360
C 11.40093 -2.61315 4.40007
C 13.71627 -0.47204 6.21416
C 15.10504 -0.35957 5.62759
N 12.87965 2.96481 2.18586
C 12.61974 3.03051 3.57230
C 12.40160 4.45267 3.90744
C 12.52798 5.11845 2.67892
C 12.78530 4.17544 1.64741
C 12.08130 4.96385 5.19671
C 12.49169 6.40342 2.01845
O 12.40022 7.57125 2.39468
C 12.45932 6.16343 0.44754
C 13.47963 7.07168 -0.18221
O 14.69205 6.91438 0.04562
O 12.87682 7.83262 -1.11385
C 13.85017 8.74263 -1.77439
H 12.17408 0.32569 -3.29524
H 13.19419 -2.97545 1.90797
H 12.46715 2.19485 5.39711
H 11.80730 5.17732 -2.45779
H 12.30971 2.83158 -3.93500
H 9.95346 3.76202 -2.34675
H 10.13479 1.99807 -2.44627
H 9.97960 2.89323 -3.92757
H 14.38286 3.79502 -3.24618
H 14.48502 5.25997 -2.20277
H 12.83081 5.88059 -4.43570
H 14.00984 4.77697 -5.04953
H 13.03329 -1.43738 -4.06601
H 12.79598 -3.06528 -3.52490
H 14.46895 -2.49240 -3.72423
H 13.80817 -4.75938 -2.48524
H 15.03082 -5.67052 -1.42616
H 15.42944 -4.13820 -2.25473
H 13.49484 -2.41215 4.37615
H 11.70463 -0.37845 5.60741
H 11.62700 -3.62900 4.77560
H 10.81663 -2.66511 3.49098
H 10.69722 -2.16886 5.10192
H 13.69104 -1.40988 6.77450
H 13.56523 0.33259 6.92601
H 15.34998 -1.39507 5.37296
H 15.80008 0.04111 6.36272
H 15.09271 0.21530 4.69431
H 12.68672 5.84898 5.33991
H 12.26531 4.22532 5.97419
H 11.01477 5.16615 5.06687
H 11.44793 6.38298 0.10337
H 13.46416 9.76032 -1.66050
H 13.96454 8.35977 -2.80594
H 14.83200 8.83737 -1.32552

