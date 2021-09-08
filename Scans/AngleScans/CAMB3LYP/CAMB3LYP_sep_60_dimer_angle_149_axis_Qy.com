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
Mg 31.53040 2.84113 2.38522
C 33.39829 0.28946 0.68845
C 33.43078 5.13132 0.71843
C 30.21824 5.09648 4.19051
C 30.59898 0.41831 4.71075
N 33.17078 2.61632 0.91864
C 33.75368 1.59881 0.30281
C 34.71666 2.00735 -0.84407
C 35.06531 3.49410 -0.34719
C 33.79797 3.78790 0.45294
C 36.25199 3.54913 0.60049
C 34.00748 2.02098 -2.24368
C 34.95364 2.10855 -3.44320
H 34.59762 1.39122 -4.66434
N 31.42758 4.80451 2.11481
C 32.26837 5.57109 1.38905
C 31.77314 6.94980 1.52348
C 30.66877 6.95637 2.52100
C 30.74423 5.59313 3.00678
C 32.28988 8.08509 0.68830
C 29.71244 8.01802 2.97325
O 28.99052 7.87146 4.02214
C 29.71395 9.27188 2.20557
N 30.72433 2.82186 4.30152
C 30.17163 3.92773 4.83019
C 29.41564 3.60727 6.09753
C 29.37255 2.01074 6.18004
C 30.22553 1.67006 4.93837
C 30.01305 4.16738 7.38818
C 27.97338 1.35074 6.31936
C 27.17561 1.55991 5.05240
N 31.75641 0.70099 2.57696
C 31.31257 -0.11147 3.64326
C 31.79911 -1.48106 3.37881
C 32.51747 -1.36482 2.17926
C 32.49823 -0.01395 1.73849
C 31.59975 -2.61167 4.22030
C 33.29664 -2.07569 1.19120
O 33.57619 -3.25699 0.99092
C 34.00825 -1.01235 0.24816
C 33.79404 -1.43956 -1.17806
O 32.65194 -1.44959 -1.67010
O 34.98807 -1.55594 -1.78725
C 34.82735 -1.96708 -3.20745
H 34.11749 5.92498 0.42279
H 29.65907 5.81387 4.79562
H 30.26822 -0.41269 5.33741
H 35.62171 1.41547 -0.90981
H 35.15451 4.18002 -1.18559
H 36.59261 2.56248 0.87732
H 35.90711 4.08888 1.47272
H 37.05235 4.15599 0.16751
H 33.46633 2.96860 -2.36166
H 33.36346 1.17012 -2.35278
H 35.99353 1.89958 -3.16407
H 34.97323 3.14291 -3.79495
H 33.21216 7.80993 0.16190
H 32.59777 8.87664 1.38004
H 31.53527 8.48460 0.00740
H 30.70979 9.70807 2.22835
H 28.90327 9.87162 2.62388
H 29.49567 9.04041 1.15706
H 28.40140 3.98413 6.04861
H 29.92405 1.63285 7.05347
H 29.30918 4.80716 7.95325
H 30.91036 4.71775 7.13824
H 30.38713 3.42032 8.08616
H 27.40830 1.82728 7.12419
H 28.01650 0.28954 6.54072
H 26.75700 2.56178 5.18694
H 26.38983 0.81208 4.96663
H 27.83181 1.59193 4.17489
H 31.33612 -3.43823 3.57393
H 30.82636 -2.42372 4.96213
H 32.59134 -2.69501 4.67315
H 35.06551 -0.99321 0.51539
H 35.42248 -2.87491 -3.34655
H 35.11012 -1.08241 -3.80840
H 33.84764 -2.30597 -3.52332

