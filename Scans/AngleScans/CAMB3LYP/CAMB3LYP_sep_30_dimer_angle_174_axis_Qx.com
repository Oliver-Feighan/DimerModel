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
Mg 15.75843 1.41648 1.20416
C 15.77932 4.12915 3.55364
C 16.40416 -0.66456 3.82846
C 16.19271 -0.92737 -0.88993
C 16.25356 3.77566 -1.31275
N 16.03378 1.80135 3.36479
C 15.95838 2.86426 4.15149
C 16.00148 2.54493 5.66996
C 16.76278 1.13283 5.59721
C 16.35092 0.70531 4.18999
C 18.27645 1.26264 5.63315
C 14.56800 2.39340 6.28980
C 14.52607 2.38828 7.81950
H 13.38697 3.01263 8.48666
N 15.74268 -0.55174 1.45798
C 16.00626 -1.22759 2.59598
C 15.88266 -2.65329 2.25501
C 15.68643 -2.76857 0.78432
C 15.89161 -1.39149 0.38227
C 15.87475 -3.73992 3.29067
C 15.37924 -3.92873 -0.11337
O 15.48739 -3.84132 -1.38753
C 15.03844 -5.19692 0.54773
N 16.37447 1.39589 -0.78139
C 16.42745 0.24638 -1.47690
C 16.62903 0.50907 -2.95025
C 16.45773 2.08864 -3.13240
C 16.26677 2.49457 -1.65451
C 17.97677 0.07393 -3.52526
C 15.40767 2.56742 -4.17203
C 14.01561 2.21513 -3.69934
N 15.79589 3.57311 1.07935
C 16.04601 4.35688 -0.06836
C 16.08072 5.76808 0.36723
C 15.86960 5.70522 1.75286
C 15.73440 4.34914 2.15567
C 16.33374 6.89112 -0.46995
C 15.74097 6.47702 2.96796
O 15.68259 7.67712 3.23272
C 15.80291 5.48318 4.20684
C 14.67453 5.83138 5.13852
O 13.49229 5.67793 4.78444
O 15.17460 6.08000 6.36265
C 14.09255 6.41936 7.32484
H 16.83013 -1.37270 4.53978
H 16.24321 -1.69012 -1.67040
H 16.30516 4.57766 -2.05221
H 16.57410 3.24581 6.26547
H 16.37491 0.43647 6.33612
H 18.59372 2.29361 5.58243
H 18.63934 0.71245 4.77487
H 18.67802 0.74968 6.51175
H 14.19597 1.38063 6.08826
H 13.90792 3.15021 5.91275
H 15.46802 2.73851 8.25884
H 14.43942 1.35386 8.16081
H 16.20296 -3.36612 4.26840
H 16.64673 -4.46171 3.00273
H 14.91414 -4.25572 3.35267
H 15.85998 -5.50034 1.19229
H 14.76441 -5.88166 -0.25720
H 14.17273 -5.03062 1.19850
H 15.87433 0.00342 -3.54003
H 17.39052 2.56333 -3.47066
H 17.88262 -0.63187 -4.37206
H 18.56243 -0.36446 -2.72817
H 18.61828 0.88580 -3.86367
H 15.55289 2.04929 -5.12313
H 15.45499 3.63214 -4.37482
H 13.90492 1.17278 -4.01293
H 13.27544 2.85235 -4.17906
H 13.95384 2.23850 -2.60513
H 15.62020 7.65483 -0.19048
H 16.24607 6.63011 -1.52255
H 17.36558 7.11643 -0.18747
H 16.77776 5.60917 4.67948
H 14.34535 7.39106 7.76024
H 14.02751 5.55884 8.01696
H 13.10690 6.61864 6.92080

