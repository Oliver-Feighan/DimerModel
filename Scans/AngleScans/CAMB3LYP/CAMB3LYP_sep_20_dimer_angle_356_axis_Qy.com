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
Mg 10.51670 0.95118 0.80579
C 10.13570 -1.28797 3.58429
C 9.11803 3.37329 2.75815
C 10.42289 2.82200 -1.75520
C 10.80483 -1.85648 -1.23864
N 9.76539 0.92432 2.88595
C 9.72884 0.03304 3.86505
C 9.30383 0.60814 5.24280
C 8.48130 1.88697 4.72618
C 9.18308 2.10436 3.38718
C 7.01578 1.59044 4.45475
C 10.53455 1.03231 6.11882
C 10.21201 1.31217 7.58831
H 11.20252 0.95213 8.59912
N 10.32922 2.92157 0.66091
C 9.75503 3.75235 1.55612
C 9.85239 5.10004 0.97438
C 10.38293 4.97596 -0.41048
C 10.37925 3.53517 -0.56608
C 9.53573 6.34582 1.74973
C 10.80941 5.98721 -1.43112
O 11.01317 5.66431 -2.65466
C 10.89093 7.38305 -0.97673
N 10.38932 0.54924 -1.23003
C 10.41970 1.54388 -2.13430
C 10.59286 1.00365 -3.53379
C 10.91641 -0.55247 -3.35664
C 10.78030 -0.66619 -1.82239
C 9.39208 1.17077 -4.46469
C 12.21491 -1.08182 -4.02475
C 13.42588 -0.49220 -3.33823
N 10.66620 -1.18312 1.10656
C 10.75286 -2.18082 0.11099
C 10.71812 -3.48730 0.79971
C 10.58947 -3.15677 2.15717
C 10.52670 -1.74450 2.30232
C 10.75289 -4.76214 0.16751
C 10.48180 -3.68172 3.49938
O 10.56248 -4.80119 4.00337
C 10.05599 -2.49887 4.47198
C 10.95338 -2.54357 5.67830
O 12.17156 -2.31862 5.56852
O 10.19661 -2.62868 6.78754
C 11.04125 -2.66474 8.01111
H 8.48538 4.14145 3.20377
H 10.50383 3.42390 -2.66330
H 10.98829 -2.77688 -1.79711
H 8.65862 -0.03839 5.82544
H 8.63137 2.74312 5.37875
H 6.79512 0.53691 4.54151
H 6.82794 1.93453 3.44611
H 6.38122 2.20018 5.10427
H 10.86941 2.02990 5.80668
H 11.31825 0.30281 6.05347
H 9.22065 0.93787 7.87124
H 10.14009 2.39354 7.72795
H 9.01296 6.11487 2.68609
H 8.80388 6.91073 1.16225
H 10.41427 6.97261 1.91694
H 9.92020 7.69785 -0.60116
H 11.29771 7.94213 -1.82162
H 11.58731 7.43602 -0.13242
H 11.42725 1.48154 -4.03230
H 10.12736 -1.18459 -3.78984
H 9.63301 1.72211 -5.39316
H 8.60392 1.67183 -3.91854
H 8.90967 0.24291 -4.76727
H 12.26228 -0.76014 -5.06803
H 12.29428 -2.16378 -4.01743
H 13.53199 0.48413 -3.82053
H 14.30285 -1.11363 -3.50788
H 13.22789 -0.31474 -2.27477
H 11.43317 -5.37631 0.74242
H 11.06856 -4.68262 -0.87066
H 9.70260 -5.05249 0.25610
H 9.00863 -2.65195 4.73498
H 10.76345 -3.56704 8.56462
H 10.87778 -1.69395 8.51567
H 12.10636 -2.81622 7.88058

