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
Mg 10.51835 0.94730 0.80586
C 9.09442 -0.97177 3.48329
C 10.47778 3.63275 2.90855
C 11.34818 2.76439 -1.65930
C 9.40131 -1.52876 -1.37749
N 9.84936 1.17553 2.90171
C 9.38306 0.36260 3.83769
C 9.29236 0.99534 5.25225
C 9.19745 2.53893 4.81958
C 9.91611 2.46030 3.47421
C 7.77326 3.00852 4.57359
C 10.57378 0.71818 6.11431
C 10.42863 1.03862 7.60361
H 11.11808 0.18724 8.56916
N 11.31506 2.76438 0.76045
C 11.21868 3.71909 1.70941
C 11.96055 4.87858 1.19040
C 12.36325 4.58809 -0.21252
C 11.65775 3.34225 -0.43669
C 12.29132 6.07635 2.03264
C 13.22848 5.31813 -1.19467
O 13.24894 5.00450 -2.43727
C 13.98005 6.47052 -0.67639
N 10.21110 0.77044 -1.24272
C 10.72238 1.67239 -2.09887
C 10.61019 1.19380 -3.52673
C 10.13423 -0.33009 -3.43293
C 9.96003 -0.44707 -1.90278
C 9.64312 1.97493 -4.41610
C 11.00998 -1.38695 -4.16009
C 12.35479 -1.49987 -3.47882
N 9.60861 -1.00274 1.00010
C 9.19792 -1.86021 -0.04401
C 8.53082 -3.02011 0.58212
C 8.57967 -2.74371 1.95679
C 9.21324 -1.48983 2.17095
C 7.93976 -4.11387 -0.11103
C 8.22985 -3.22255 3.27476
O 7.75467 -4.26546 3.72229
C 8.43465 -2.03735 4.31386
C 9.19650 -2.57918 5.49222
O 10.36982 -2.96990 5.36076
O 8.49429 -2.34587 6.61598
C 9.21429 -2.85543 7.81337
H 10.29981 4.58595 3.40717
H 11.71218 3.29954 -2.53941
H 9.11286 -2.38990 -1.98405
H 8.41388 0.71374 5.82036
H 9.74582 3.17656 5.50810
H 7.06707 2.19263 4.61574
H 7.77690 3.45521 3.58795
H 7.51642 3.81336 5.26828
H 11.35240 1.44206 5.84138
H 10.90250 -0.29568 5.99319
H 9.38058 1.17925 7.89474
H 10.89292 2.00876 7.79668
H 11.72233 6.07810 2.97053
H 11.92760 6.95725 1.49264
H 13.36395 6.18608 2.20606
H 13.28591 7.19684 -0.26035
H 14.60770 6.80629 -1.50415
H 14.61397 6.13152 0.15056
H 11.57166 1.23169 -4.02397
H 9.13712 -0.47344 -3.87459
H 10.12217 2.38925 -5.32325
H 9.19918 2.76544 -3.82572
H 8.76961 1.41734 -4.74970
H 11.20808 -1.07237 -5.18769
H 10.55192 -2.36929 -4.20666
H 12.92330 -0.67386 -3.91653
H 12.81762 -2.45918 -3.70140
H 12.26847 -1.30708 -2.40316
H 8.23443 -5.01192 0.41544
H 8.25410 -4.14127 -1.15228
H 6.88116 -3.86070 -0.00835
H 7.44553 -1.67549 4.59718
H 8.53196 -3.53724 8.33033
H 9.54477 -1.95713 8.36799
H 10.07047 -3.49872 7.64730

