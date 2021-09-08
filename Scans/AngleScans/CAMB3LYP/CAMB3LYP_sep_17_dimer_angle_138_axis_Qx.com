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
Mg 8.92571 0.80639 0.68463
C 10.36097 2.91013 3.21314
C 8.03350 -1.33577 3.18458
C 8.05515 -1.27602 -1.54542
C 10.89341 2.49702 -1.64252
N 9.21493 0.88727 2.87600
C 9.71971 1.76504 3.72985
C 9.45574 1.43454 5.22338
C 9.24704 -0.15174 5.08540
C 8.76739 -0.21005 3.63654
C 10.54152 -0.94058 5.19209
C 8.16508 2.13811 5.77219
C 8.01597 2.11098 7.29494
H 7.41430 3.26571 7.95606
N 7.74008 -0.78075 0.80229
C 7.47281 -1.51712 1.90120
C 6.56202 -2.58426 1.45868
C 6.44408 -2.51573 -0.02325
C 7.44670 -1.51275 -0.32137
C 5.84235 -3.48890 2.41635
C 5.58179 -3.24295 -1.01006
O 5.81382 -3.19642 -2.26975
C 4.51461 -4.08616 -0.45182
N 9.55630 0.48848 -1.27077
C 8.97595 -0.44881 -2.04058
C 9.40076 -0.30970 -3.48300
C 10.20231 1.07184 -3.56347
C 10.17785 1.46594 -2.07033
C 10.27492 -1.43760 -4.03077
C 9.71238 2.10980 -4.61002
C 8.34817 2.63152 -4.21995
N 10.22986 2.52878 0.70935
C 10.97551 3.04967 -0.37082
C 11.79915 4.15456 0.16140
C 11.49028 4.18509 1.52966
C 10.55637 3.15777 1.83284
C 12.72332 4.93788 -0.58590
C 11.74996 4.84594 2.78845
O 12.38724 5.84068 3.13211
C 11.12614 3.96866 3.95790
C 10.35175 4.88594 4.86423
O 9.33410 5.47020 4.45237
O 10.81109 4.75356 6.12195
C 10.06664 5.63544 7.06000
H 7.90960 -2.18064 3.86261
H 7.70585 -1.89705 -2.37354
H 11.45965 3.13694 -2.32262
H 10.28497 1.64401 5.88842
H 8.47154 -0.50812 5.75842
H 11.40580 -0.29399 5.22509
H 10.57458 -1.57182 4.31375
H 10.50013 -1.61877 6.04923
H 7.28582 1.54636 5.48673
H 8.10414 3.14987 5.42110
H 8.94893 1.82456 7.79543
H 7.31437 1.31657 7.56082
H 6.25454 -3.41124 3.42999
H 6.06290 -4.51764 2.11153
H 4.76043 -3.34074 2.40351
H 4.95205 -4.83555 0.20356
H 3.95108 -4.45213 -1.31227
H 3.86605 -3.46180 0.17280
H 8.53873 -0.25441 -4.13641
H 11.25795 0.91551 -3.82997
H 9.84721 -1.92531 -4.92703
H 10.43173 -2.16159 -3.24224
H 11.29336 -1.15004 -4.28626
H 9.59547 1.63568 -5.58761
H 10.38984 2.94755 -4.73736
H 7.67061 1.86530 -4.60842
H 8.16004 3.59719 -4.68485
H 8.23174 2.65267 -3.13016
H 12.57512 5.96625 -0.28426
H 12.57679 4.81175 -1.65660
H 13.66705 4.50252 -0.24665
H 11.95168 3.48082 4.47757
H 10.80847 6.25704 7.57087
H 9.45872 4.95769 7.68835
H 9.41810 6.39002 6.63053

