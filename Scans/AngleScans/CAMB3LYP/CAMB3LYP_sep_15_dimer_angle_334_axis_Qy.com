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
Mg 7.88564 0.71454 0.60793
C 6.61473 -1.64014 2.99944
C 5.75616 3.01049 1.96016
C 8.65273 2.67647 -1.76479
C 9.03539 -2.00193 -1.24804
N 6.42448 0.58264 2.26318
C 6.07214 -0.34578 3.13967
C 5.14421 0.16152 4.27596
C 4.51057 1.42935 3.52099
C 5.64444 1.72173 2.54042
C 3.26291 1.09208 2.72165
C 5.94561 0.59514 5.55335
C 5.09288 0.80771 6.80611
H 5.65813 0.44464 8.10270
N 7.67304 2.68121 0.44801
C 6.77200 3.45686 1.08657
C 7.01352 4.82819 0.61190
C 8.02110 4.77529 -0.48204
C 8.14193 3.34218 -0.66005
C 6.37642 6.03187 1.24332
C 8.74530 5.83876 -1.25064
O 9.39930 5.56974 -2.31963
C 8.58903 7.21841 -0.76728
N 8.53431 0.38570 -1.33975
C 8.84860 1.41429 -2.14669
C 9.54898 0.93381 -3.39529
C 9.85682 -0.61559 -3.14598
C 9.17177 -0.79199 -1.77303
C 8.76803 1.09382 -4.69947
C 11.33312 -1.07341 -3.29966
C 12.17813 -0.46808 -2.20197
N 8.01330 -1.42306 0.89505
C 8.50616 -2.37872 -0.02031
C 8.28156 -3.71074 0.57779
C 7.64774 -3.43652 1.79908
C 7.47035 -2.03384 1.94225
C 8.60557 -4.95875 -0.02522
C 7.07883 -4.01534 2.99491
O 7.02063 -5.14957 3.46808
C 6.27078 -2.88581 3.76796
C 6.66301 -2.94462 5.21890
O 7.82437 -2.67308 5.57119
O 5.55634 -3.09813 5.96863
C 5.89282 -3.15082 7.41633
H 4.96899 3.73857 2.15800
H 9.03367 3.31478 -2.56538
H 9.45392 -2.89353 -1.71980
H 4.36086 -0.52883 4.56501
H 4.37025 2.26492 4.20179
H 3.07510 0.02893 2.69743
H 3.44315 1.46728 1.72276
H 2.40639 1.65418 3.10473
H 6.32494 1.61496 5.40901
H 6.73152 -0.10343 5.76544
H 4.08550 0.38843 6.69505
H 4.92442 1.87983 6.93342
H 5.55744 5.74746 1.91546
H 5.88629 6.59259 0.44006
H 7.10181 6.68231 1.73667
H 7.53469 7.48454 -0.76932
H 9.25138 7.82298 -1.38990
H 8.92305 7.26379 0.27531
H 10.48497 1.45919 -3.54014
H 9.31257 -1.25816 -3.85360
H 9.30740 1.68801 -5.46110
H 7.81189 1.54590 -4.47152
H 8.47435 0.16173 -5.17922
H 11.74561 -0.71102 -4.24448
H 11.45452 -2.15143 -3.28767
H 12.40852 0.52892 -2.58927
H 13.08398 -1.05151 -2.04999
H 11.59510 -0.33805 -1.28284
H 9.05471 -5.56999 0.74624
H 9.27666 -4.82897 -0.87172
H 7.61097 -5.28890 -0.33674
H 5.20841 -3.08537 3.62265
H 5.47335 -4.08261 7.80813
H 5.51035 -2.20614 7.84650
H 6.93716 -3.25980 7.68451

