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
Mg 15.77587 1.42363 1.19609
C 17.80026 -0.69789 3.26493
C 16.82008 4.04163 3.11711
C 13.66327 3.24608 -0.31489
C 14.04345 -1.43245 0.20250
N 17.13574 1.53147 2.93718
C 17.89524 0.66973 3.59665
C 18.86729 1.31931 4.61778
C 18.06494 2.68588 4.87779
C 17.31638 2.78618 3.55029
C 17.04169 2.57949 5.99627
C 20.28458 1.59346 4.00285
C 21.37221 1.93854 5.02252
H 22.73266 1.46919 4.77481
N 15.66312 3.39879 1.03922
C 16.16200 4.31304 1.89744
C 15.79621 5.62607 1.34398
C 14.89911 5.40957 0.17648
C 14.68565 3.97925 0.26972
C 16.35328 6.91582 1.87276
C 14.31790 6.33837 -0.84614
O 13.37295 5.97046 -1.63007
C 14.82292 7.71909 -0.83994
N 13.96376 1.00678 0.26613
C 13.27001 1.97307 -0.36071
C 12.14795 1.39071 -1.18668
C 12.38312 -0.19111 -1.17635
C 13.60259 -0.25948 -0.23103
C 10.73289 1.69779 -0.69709
C 12.47718 -0.89767 -2.55642
C 13.73284 -0.46022 -3.27568
N 15.99127 -0.70632 1.48754
C 15.13968 -1.72381 1.00424
C 15.63147 -3.00236 1.55732
C 16.72996 -2.63441 2.34879
C 16.89777 -1.22383 2.30907
C 15.04441 -4.28148 1.34441
C 17.77925 -3.11745 3.21733
O 18.18634 -4.22896 3.55282
C 18.44182 -1.87236 3.95022
C 19.93577 -2.01354 3.84682
O 20.50139 -1.95224 2.74096
O 20.46813 -1.97948 5.08205
C 21.94918 -2.10489 5.03100
H 16.90415 4.89371 3.79222
H 12.97133 3.81640 -0.93893
H 13.61727 -2.37846 -0.13823
H 18.98219 0.77316 5.54638
H 18.74405 3.52608 4.99702
H 16.93965 1.56567 6.35389
H 16.10785 2.92778 5.57474
H 17.28846 3.27833 6.80067
H 20.25340 2.53296 3.43618
H 20.60446 0.76644 3.39923
H 21.06509 1.70244 6.04862
H 21.51216 3.02224 5.02727
H 16.85539 6.77165 2.83739
H 15.49949 7.56162 2.10482
H 16.99539 7.42476 1.15072
H 14.64331 8.16494 0.13547
H 14.35487 8.20531 -1.69808
H 15.90944 7.69483 -0.97963
H 12.19519 1.74635 -2.20863
H 11.56316 -0.72166 -0.67037
H 10.10575 2.19665 -1.45992
H 10.80473 2.30726 0.19393
H 10.16963 0.83614 -0.34282
H 11.63754 -0.60304 -3.19073
H 12.46548 -1.98040 -2.48903
H 13.43523 0.48528 -3.73887
H 14.01995 -1.19415 -4.02599
H 14.53737 -0.23999 -2.56442
H 15.85748 -4.96955 1.15478
H 14.33794 -4.26198 0.51719
H 14.54627 -4.42998 2.30622
H 18.10040 -1.88217 4.98604
H 22.22009 -2.95326 5.66704
H 22.34316 -1.11247 5.32019
H 22.39520 -2.39700 4.08750

