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
Mg 13.14662 1.17633 0.99578
C 13.19266 1.38221 -2.58676
C 12.59128 -2.17101 0.64722
C 12.65271 1.03306 4.12674
C 12.56430 4.61202 1.04712
N 12.93927 -0.11098 -0.79086
C 13.03657 0.06775 -2.09965
C 13.04295 -1.24423 -2.92926
C 12.28395 -2.19517 -1.88116
C 12.65194 -1.47417 -0.58601
C 10.77182 -2.16830 -2.02973
C 14.49594 -1.75827 -3.22354
C 14.58672 -2.85721 -4.28460
H 15.74467 -2.87184 -5.17421
N 13.17638 -0.37680 2.23087
C 12.95130 -1.67004 1.91755
C 13.06822 -2.41610 3.18004
C 13.21771 -1.43740 4.29114
C 12.99567 -0.19470 3.57946
C 13.11246 -3.91546 3.23807
C 13.49956 -1.59480 5.75469
O 13.35051 -0.62332 6.57740
C 13.86510 -2.94390 6.21055
N 12.46755 2.56974 2.38128
C 12.39583 2.26586 3.68912
C 12.14650 3.49992 4.52299
C 12.30717 4.73552 3.52047
C 12.54401 3.96389 2.20370
C 10.78240 3.57513 5.20851
C 13.32205 5.84057 3.92242
C 14.72955 5.29106 3.87379
N 13.09874 2.76777 -0.46473
C 12.80974 4.13033 -0.23238
C 12.78473 4.80068 -1.54860
C 13.04018 3.76900 -2.46450
C 13.19224 2.53864 -1.76981
C 12.50174 6.17704 -1.77593
C 13.20525 3.43914 -3.86187
O 13.26846 4.08715 -4.90581
C 13.18588 1.85703 -4.01333
C 14.34240 1.46010 -4.88937
O 15.51318 1.63650 -4.50889
O 13.88093 0.74346 -5.93045
C 14.99214 0.31735 -6.82233
H 12.19037 -3.18498 0.65101
H 12.57960 1.05966 5.21644
H 12.48671 5.69964 0.98600
H 12.48755 -1.19695 -3.85827
H 12.69730 -3.20038 -1.88866
H 10.45000 -1.42141 -2.74017
H 10.38335 -1.94561 -1.04457
H 10.40005 -3.16555 -2.28195
H 14.86437 -2.31030 -2.34932
H 15.14137 -0.94412 -3.49038
H 13.65824 -2.95158 -4.86077
H 14.68732 -3.82054 -3.77874
H 12.81452 -4.36397 2.28232
H 12.33384 -4.23137 3.94079
H 14.07610 -4.29534 3.58413
H 13.06547 -3.63787 5.96264
H 14.11533 -2.83729 7.26786
H 14.75065 -3.27281 5.65537
H 12.88350 3.58913 5.31165
H 11.36264 5.28545 3.39691
H 10.85157 3.69258 6.30647
H 10.22381 2.68368 4.95589
H 10.12798 4.36744 4.84908
H 13.14808 6.15758 4.95363
H 13.26509 6.72672 3.29913
H 14.83328 4.79220 4.84202
H 15.45210 6.09747 3.76585
H 14.82616 4.52461 3.09606
H 13.22155 6.52677 -2.50399
H 12.55653 6.75180 -0.85375
H 11.47878 6.10580 -2.15505
H 12.22626 1.58170 -4.45253
H 14.75047 0.67611 -7.82760
H 15.08183 -0.77678 -6.68600
H 15.96378 0.77044 -6.66398

