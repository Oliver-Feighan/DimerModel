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
Mg 15.76412 1.42173 1.18755
C 17.63201 -1.12994 -0.50922
C 17.66450 3.71192 -0.47923
C 14.45196 3.67707 2.99285
C 14.83270 -1.00110 3.51308
N 17.40450 1.19692 -0.27903
C 17.98740 0.17941 -0.89486
C 18.95038 0.58795 -2.04174
C 19.29903 2.07470 -1.54485
C 18.03169 2.36850 -0.74472
C 20.48571 2.12973 -0.59717
C 18.24120 0.60158 -3.44134
C 19.18736 0.68914 -4.64087
H 18.83134 -0.02818 -5.86200
N 15.66130 3.38511 0.91714
C 16.50209 4.15169 0.19139
C 16.00686 5.53040 0.32582
C 14.90249 5.53697 1.32334
C 14.97795 4.17373 1.80911
C 16.52360 6.66569 -0.50937
C 13.94616 6.59862 1.77559
O 13.22424 6.45206 2.82448
C 13.94767 7.85248 1.00791
N 14.95806 1.40246 3.10386
C 14.40536 2.50833 3.63252
C 13.64936 2.18787 4.89987
C 13.60627 0.59134 4.98238
C 14.45925 0.25066 3.74071
C 14.24677 2.74798 6.19052
C 12.20710 -0.06866 5.12170
C 11.40933 0.14051 3.85474
N 15.99014 -0.71841 1.37930
C 15.54629 -1.53087 2.44560
C 16.03283 -2.90046 2.18115
C 16.75119 -2.78422 0.98159
C 16.73195 -1.43335 0.54083
C 15.83348 -4.03107 3.02263
C 17.53036 -3.49509 -0.00646
O 17.80991 -4.67640 -0.20674
C 18.24197 -2.43175 -0.94950
C 18.02776 -2.85896 -2.37573
O 16.88566 -2.86899 -2.86777
O 19.22180 -2.97534 -2.98492
C 19.06107 -3.38648 -4.40511
H 18.35122 4.50558 -0.77488
H 13.89280 4.39447 3.59796
H 14.50194 -1.83210 4.13975
H 19.85543 -0.00393 -2.10748
H 19.38823 2.76062 -2.38325
H 20.82633 1.14308 -0.32035
H 20.14083 2.66948 0.27505
H 21.28608 2.73659 -1.03015
H 17.70005 1.54920 -3.55933
H 17.59718 -0.24928 -3.55044
H 20.22725 0.48018 -4.36173
H 19.20695 1.72351 -4.99261
H 17.44588 6.39053 -1.03576
H 16.83149 7.45724 0.18238
H 15.76899 7.06520 -1.19027
H 14.94351 8.28867 1.03069
H 13.13699 8.45222 1.42622
H 13.72939 7.62101 -0.04061
H 12.63512 2.56473 4.85094
H 14.15777 0.21345 5.85581
H 13.54290 3.38776 6.75559
H 15.14408 3.29835 5.94057
H 14.62085 2.00092 6.88849
H 11.64202 0.40788 5.92652
H 12.25022 -1.12986 5.34305
H 10.99073 1.14238 3.98927
H 10.62356 -0.60732 3.76897
H 12.06553 0.17253 2.97723
H 15.56985 -4.85763 2.37627
H 15.06008 -3.84312 3.76446
H 16.82506 -4.11441 3.47549
H 19.29923 -2.41261 -0.68227
H 19.65620 -4.29431 -4.54422
H 19.34384 -2.50181 -5.00606
H 18.08136 -3.72537 -4.72098

