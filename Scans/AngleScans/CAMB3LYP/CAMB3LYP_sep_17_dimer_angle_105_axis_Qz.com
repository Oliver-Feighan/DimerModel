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
Mg 8.94190 0.79766 0.68167
C 9.14514 -1.33818 -2.19513
C 8.17039 -1.98048 2.50412
C 8.31458 2.70499 3.13850
C 8.59737 3.48497 -1.51043
N 8.71740 -1.32851 0.11485
C 8.87905 -2.03152 -0.99600
C 8.82767 -3.57009 -0.79667
C 7.96336 -3.60856 0.55624
C 8.32838 -2.23626 1.11858
C 6.46387 -3.64109 0.31171
C 10.25215 -4.19343 -0.58619
C 10.30890 -5.71887 -0.69381
H 11.49767 -6.33150 -1.28031
N 8.81464 0.39455 2.62072
C 8.51322 -0.79281 3.18693
C 8.52766 -0.56102 4.63959
C 8.70009 0.89743 4.88020
C 8.59326 1.40243 3.52608
C 8.46558 -1.67650 5.64223
C 8.91167 1.70535 6.12466
O 8.79753 2.98199 6.12618
C 9.16438 0.95096 7.36100
N 8.30617 2.77315 0.80662
C 8.16139 3.37874 1.99836
C 7.96504 4.86694 1.83371
C 8.25083 5.16976 0.28963
C 8.48606 3.72798 -0.21170
C 6.58313 5.40121 2.20956
C 9.32294 6.24744 -0.02951
C 10.68982 5.75533 0.38867
N 9.06274 1.08642 -1.45610
C 8.85977 2.28945 -2.16715
C 8.93384 1.96260 -3.60593
C 9.15360 0.57715 -3.63264
C 9.19220 0.07269 -2.30466
C 8.75618 2.88189 -4.67807
C 9.35110 -0.57499 -4.48251
O 9.50069 -0.74719 -5.69152
C 9.22825 -1.88631 -3.59272
C 10.38856 -2.78379 -3.92546
O 11.55262 -2.43730 -3.65802
O 9.92052 -3.98756 -4.30263
C 11.03428 -4.91635 -4.63234
H 7.70038 -2.74559 3.12266
H 8.20002 3.42479 3.95228
H 8.59780 4.28321 -2.25575
H 8.31442 -4.11291 -1.58134
H 8.30606 -4.39673 1.22160
H 6.22319 -3.51328 -0.73322
H 6.05272 -2.82886 0.89676
H 6.03428 -4.55799 0.72533
H 10.54625 -4.06809 0.46383
H 10.96257 -3.75689 -1.26099
H 9.39987 -6.13432 -1.14549
H 10.32225 -6.13769 0.31536
H 8.17571 -2.62430 5.17207
H 7.63951 -1.44780 6.32426
H 9.38600 -1.77347 6.22197
H 8.32915 0.28156 7.55270
H 9.37896 1.70251 8.12323
H 10.04649 0.31855 7.21083
H 8.67436 5.41972 2.43737
H 7.35234 5.53917 -0.22623
H 6.61636 6.19212 2.98250
H 5.97468 4.57046 2.54149
H 6.00005 5.79775 1.38012
H 9.13030 7.15562 0.54690
H 9.35245 6.53063 -1.07637
H 10.71999 5.98899 1.45712
H 11.47023 6.28548 -0.15349
H 10.76409 4.66637 0.28717
H 9.52698 2.66457 -5.40548
H 8.81376 3.91188 -4.33230
H 7.74669 2.61292 -5.00035
H 8.27017 -2.35220 -3.82639
H 10.85838 -5.27747 -5.65035
H 11.04245 -5.67178 -3.82426
H 12.02787 -4.49400 -4.72648

