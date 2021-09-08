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
Mg 15.77565 1.42348 1.19486
C 18.01534 -0.70780 3.01708
C 17.08727 4.04387 2.94064
C 13.51450 3.26224 -0.05957
C 13.89454 -1.41628 0.45803
N 17.34455 1.52782 2.75046
C 18.16931 0.66133 3.31907
C 19.27034 1.30624 4.20303
C 18.52514 2.68670 4.54636
C 17.61732 2.78545 3.32212
C 17.64901 2.60568 5.78538
C 20.60278 1.55397 3.41230
C 21.81424 1.89272 4.28371
H 23.12651 1.40133 3.87266
N 15.67030 3.39851 1.03142
C 16.28501 4.31340 1.81013
C 15.87009 5.62637 1.29233
C 14.83080 5.41203 0.24895
C 14.61178 3.98590 0.38410
C 16.50617 6.91278 1.73276
C 14.13822 6.33957 -0.70300
O 13.09759 5.97812 -1.35818
C 14.65829 7.71280 -0.77549
N 13.85580 1.02434 0.50404
C 13.10175 1.99471 -0.04159
C 11.87730 1.42103 -0.71389
C 12.09092 -0.16388 -0.71555
C 13.41835 -0.24113 0.07017
C 10.53905 1.75323 -0.05425
C 12.00169 -0.88462 -2.08859
C 13.16286 -0.47223 -2.96438
N 15.99767 -0.70656 1.48065
C 15.07875 -1.71604 1.11924
C 15.61906 -2.99637 1.62048
C 16.81296 -2.63705 2.26389
C 16.99315 -1.22949 2.18779
C 14.99302 -4.26875 1.49706
C 17.95645 -3.12715 2.99938
O 18.38765 -4.24128 3.29353
C 18.72219 -1.88502 3.62956
C 20.18936 -2.04888 3.34130
O 20.61250 -2.00617 2.17266
O 20.87292 -2.01101 4.49962
C 22.33405 -2.15841 4.26473
H 17.26668 4.90092 3.59038
H 12.75734 3.83670 -0.59826
H 13.41647 -2.35914 0.18393
H 19.49361 0.76718 5.11591
H 19.22493 3.51802 4.57019
H 17.57922 1.59684 6.16421
H 16.67435 2.96355 5.48038
H 18.00401 3.30835 6.54470
H 20.51319 2.48848 2.84361
H 20.83339 0.71677 2.78257
H 21.63519 1.67073 5.34277
H 21.96803 2.97427 4.25885
H 17.12339 6.77035 2.62839
H 15.69687 7.57307 2.06281
H 17.05929 7.40554 0.93032
H 14.60842 8.17032 0.20972
H 14.09275 8.19772 -1.57355
H 15.71827 7.67143 -1.04998
H 11.80064 1.76636 -1.73760
H 11.33396 -0.67769 -0.10492
H 9.82782 2.25400 -0.73795
H 10.73020 2.36990 0.81392
H 10.01333 0.90321 0.37738
H 11.09306 -0.58376 -2.61589
H 11.98418 -1.96639 -2.00824
H 12.82204 0.47313 -3.39708
H 13.34379 -1.21725 -3.73655
H 14.05312 -0.25708 -2.36206
H 15.76668 -4.97032 1.21467
H 14.18866 -4.24672 0.76474
H 14.61758 -4.40098 2.51531
H 18.51335 -1.88017 4.70005
H 22.67137 -3.00465 4.87118
H 22.77434 -1.16917 4.49122
H 22.65425 -2.46579 3.27607

