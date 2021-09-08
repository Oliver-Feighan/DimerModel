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
Mg 7.87553 0.70633 0.60540
C 7.73967 3.42787 2.94073
C 8.66652 -1.31547 3.23623
C 8.44396 -1.61709 -1.47932
C 8.20626 3.07804 -1.92666
N 8.13897 1.11968 2.76225
C 8.00152 2.18006 3.54393
C 8.07365 1.87247 5.06372
C 8.92185 0.51066 4.99344
C 8.52934 0.05033 3.59105
C 10.42454 0.73560 5.01920
C 6.65627 1.63452 5.69333
C 6.62384 1.63522 7.22327
H 5.45170 2.19037 7.89436
N 7.98516 -1.25756 0.86941
C 8.29751 -1.90920 2.00919
C 8.26182 -3.34171 1.67633
C 8.06450 -3.47722 0.20751
C 8.18023 -2.09220 -0.20289
C 8.32846 -4.42095 2.71759
C 7.82558 -4.65935 -0.68225
O 7.92043 -4.57234 -1.95750
C 7.56918 -5.94278 -0.01252
N 8.47982 0.71356 -1.38386
C 8.60089 -0.43417 -2.07378
C 8.77677 -0.16746 -3.54969
C 8.50534 1.39718 -3.73887
C 8.29802 1.79845 -2.26191
C 10.14578 -0.52016 -4.13092
C 7.42108 1.80323 -4.77432
C 6.05677 1.36671 -4.29109
N 7.77648 2.86033 0.46927
C 7.96996 3.65192 -0.68400
C 7.91840 5.06489 -0.25588
C 7.71990 4.99654 1.13135
C 7.67269 3.63688 1.54196
C 8.09527 6.19698 -1.10040
C 7.55019 5.76542 2.34326
O 7.41799 6.96091 2.60221
C 7.68191 4.78429 3.58681
C 6.53942 5.06598 4.52376
O 5.36709 4.83653 4.17791
O 7.03013 5.35232 5.74342
C 5.93461 5.62827 6.71064
H 9.14042 -1.99148 3.94849
H 8.53770 -2.37946 -2.25616
H 8.20289 3.87760 -2.67054
H 8.60457 2.61125 5.65202
H 8.58297 -0.20463 5.73835
H 10.67600 1.78419 4.96119
H 10.81621 0.20459 4.16149
H 10.86281 0.25377 5.89789
H 6.34751 0.59926 5.49933
H 5.94764 2.34622 5.31655
H 7.54448 2.04641 7.65487
H 6.60448 0.59930 7.57042
H 8.63831 -4.02186 3.69131
H 9.14260 -5.09435 2.42851
H 7.40260 -4.99579 2.78827
H 8.41201 -6.19038 0.62841
H 7.33400 -6.64784 -0.81219
H 6.69861 -5.82766 0.64281
H 8.05188 -0.72282 -4.13211
H 9.40439 1.92771 -4.08542
H 10.09119 -1.23514 -4.97347
H 10.76259 -0.91644 -3.33528
H 10.73291 0.32858 -4.47751
H 7.59295 1.29001 -5.72365
H 7.40010 2.86768 -4.98288
H 6.01003 0.31774 -4.59862
H 5.27514 1.95346 -4.76943
H 6.00017 1.39219 -3.19665
H 7.33677 6.91584 -0.82037
H 8.01794 5.92517 -2.15107
H 9.11255 6.48830 -0.82557
H 8.64969 4.97394 4.05265
H 6.12836 6.61633 7.13945
H 5.92796 4.76919 7.40756
H 4.93599 5.76293 6.31179

