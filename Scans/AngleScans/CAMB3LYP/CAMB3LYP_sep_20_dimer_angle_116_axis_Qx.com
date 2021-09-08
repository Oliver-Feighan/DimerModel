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
Mg 10.50345 0.95118 0.80422
C 12.50515 2.31411 3.45274
C 8.78180 -0.77050 3.19301
C 9.02167 -0.60768 -1.52852
C 13.05069 1.85210 -1.39705
N 10.71025 0.86612 3.00473
C 11.46647 1.47443 3.90609
C 11.03720 1.22843 5.37740
C 10.26532 -0.16486 5.17294
C 9.85890 -0.00663 3.70923
C 11.17295 -1.37824 5.28814
C 10.07432 2.34520 5.91380
C 9.86217 2.33793 7.42928
H 9.70063 3.61672 8.11563
N 8.81312 -0.08847 0.82588
C 8.24791 -0.70070 1.88748
C 7.02766 -1.34507 1.37763
C 7.00529 -1.20173 -0.10347
C 8.31819 -0.63282 -0.33310
C 5.98621 -1.94316 2.27824
C 5.97795 -1.53494 -1.14257
O 6.26327 -1.54667 -2.39209
C 4.65301 -1.93793 -0.64897
N 11.05395 0.47069 -1.14217
C 10.20211 -0.16706 -1.96410
C 10.70821 -0.15949 -3.38684
C 11.96456 0.83004 -3.39067
C 12.02445 1.16881 -1.88500
C 11.12799 -1.51673 -3.95081
C 11.93523 2.00073 -4.41110
C 10.84376 2.97956 -4.04238
N 12.34770 2.06926 0.93345
C 13.27718 2.30422 -0.10335
C 14.42662 3.01365 0.49492
C 14.09375 3.12270 1.85352
C 12.83550 2.50568 2.08932
C 15.60458 3.41841 -0.19407
C 14.52554 3.61009 3.14375
O 15.46922 4.29059 3.54393
C 13.57432 2.99687 4.25979
C 13.15471 4.11285 5.17686
O 12.44169 5.04140 4.75733
O 13.48007 3.78963 6.44179
C 13.07392 4.86090 7.39007
H 8.32743 -1.52620 3.83438
H 8.50328 -1.03550 -2.38976
H 13.84061 2.25414 -2.03507
H 11.85673 1.10068 6.07441
H 9.38568 -0.22605 5.80843
H 12.21233 -1.09741 5.37153
H 11.00806 -1.95547 4.38777
H 10.84907 -2.01390 6.11721
H 9.05176 2.12701 5.57988
H 10.40475 3.31621 5.59983
H 10.60244 1.71519 7.94614
H 8.90699 1.85243 7.64364
H 6.35528 -2.04789 3.30596
H 5.82520 -2.97289 1.94119
H 5.03630 -1.40577 2.23874
H 4.75614 -2.81161 -0.00983
H 4.03085 -2.04890 -1.53919
H 4.25425 -1.13366 -0.02068
H 9.95516 0.22605 -4.06313
H 12.89872 0.30157 -3.63127
H 10.58869 -1.79008 -4.87737
H 10.97424 -2.26646 -3.18608
H 12.19052 -1.61935 -4.16463
H 11.69302 1.62724 -5.40910
H 12.87815 2.53202 -4.48615
H 9.94864 2.52734 -4.47980
H 11.04384 3.95749 -4.47558
H 10.69782 3.01563 -2.95653
H 15.83278 4.42116 0.14195
H 15.46683 3.38143 -1.27269
H 16.30678 2.65726 0.15640
H 14.13989 2.22620 4.78494
H 13.97039 5.15206 7.94618
H 12.23353 4.44037 7.97373
H 12.76717 5.81184 6.97047

