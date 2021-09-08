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
Mg 31.53992 2.84330 2.38917
C 30.96707 6.31811 3.07956
C 30.84300 3.53127 -0.87818
C 31.59495 -0.27829 1.82344
C 31.06757 2.23310 5.78774
N 30.99033 4.72078 1.35705
C 30.86305 5.99061 1.71159
C 30.68981 6.97559 0.52442
C 30.09968 5.94635 -0.55784
C 30.70668 4.65194 -0.02050
C 28.58605 5.82137 -0.50664
C 32.05509 7.58695 0.05088
C 31.93033 8.78417 -0.89403
H 32.91492 9.85770 -0.79107
N 31.74266 1.85079 0.68260
C 31.43195 2.28897 -0.55543
C 31.74686 1.17759 -1.46638
C 32.11364 -0.01010 -0.64809
C 31.80443 0.47941 0.68039
C 31.76039 1.32847 -2.95986
C 32.64165 -1.36693 -1.00362
O 32.66415 -2.32497 -0.15245
C 33.04473 -1.56749 -2.40319
N 31.14994 1.18683 3.58347
C 31.29889 -0.06089 3.10500
C 31.23165 -1.08244 4.21509
C 31.24834 -0.24093 5.57493
C 31.23160 1.18781 4.98859
C 30.01022 -2.00130 4.20048
C 32.34690 -0.60141 6.61219
C 33.70922 -0.23293 6.07027
N 31.28111 4.04738 4.16435
C 31.07354 3.58840 5.48357
C 30.83772 4.77845 6.32654
C 30.90213 5.85413 5.42807
C 31.13983 5.36723 4.11432
C 30.55641 4.76703 7.72184
C 30.81253 7.29002 5.29056
O 30.71035 8.23302 6.07417
C 30.72476 7.64505 3.74375
C 31.70007 8.75599 3.46610
O 32.92350 8.56760 3.58579
O 31.04540 9.78556 2.89885
C 31.97258 10.90787 2.59458
H 30.42189 3.59201 -1.88211
H 31.71339 -1.35750 1.70201
H 31.00965 2.13420 6.87380
H 29.98262 7.77757 0.69906
H 30.47837 6.16033 -1.55391
H 28.16557 6.35814 0.33075
H 28.38111 4.76245 -0.41867
H 28.14949 6.13471 -1.45930
H 32.55518 6.87271 -0.61583
H 32.66576 7.85245 0.89184
H 30.91354 9.19518 -0.90657
H 32.09165 8.43682 -1.91742
H 31.28873 2.26791 -3.27354
H 31.10783 0.54805 -3.36601
H 32.75905 1.21125 -3.38602
H 32.19588 -1.37409 -3.05473
H 33.47796 -2.56883 -2.44190
H 33.81099 -0.82626 -2.65623
H 32.09680 -1.73373 4.19422
H 30.31178 -0.36000 6.13929
H 30.27259 -3.07561 4.16752
H 29.39253 -1.73455 3.35322
H 29.32466 -1.87150 5.03616
H 32.36365 -1.68003 6.78671
H 32.20604 -0.12155 7.57494
H 33.96655 -1.09255 5.44424
H 34.42327 -0.10780 6.88177
H 33.64842 0.64470 5.41628
H 31.14758 5.55731 8.16511
H 30.78626 3.80071 8.16586
H 29.48153 4.96613 7.70572
H 29.69631 7.94108 3.53342
H 31.56928 11.79799 3.08724
H 32.05532 10.93642 1.49181
H 32.96871 10.86420 3.01914

