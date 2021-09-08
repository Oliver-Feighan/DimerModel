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
Mg 31.53886 2.84579 2.40021
C 31.16351 1.82395 5.81986
C 30.31101 5.88894 3.33094
C 31.43252 3.63162 -0.67202
C 31.63862 -0.53178 1.54710
N 30.86066 3.63350 4.35260
C 30.82173 3.17364 5.59411
C 30.46955 4.24110 6.66465
C 29.68294 5.27643 5.72219
C 30.34561 4.94360 4.38699
C 28.19778 4.97605 5.60633
C 31.74678 4.89681 7.29769
C 31.48794 5.71951 8.56177
H 32.49763 5.71092 9.61672
N 31.42803 4.62766 1.53347
C 30.92049 5.76023 2.06354
C 31.05298 6.78779 1.01917
C 31.52894 6.13025 -0.22816
C 31.46022 4.73725 0.16546
C 30.81555 8.24714 1.27871
C 31.96072 6.66577 -1.55969
O 32.10782 5.90080 -2.57754
C 32.11580 8.12427 -1.66072
N 31.32349 1.72240 0.66413
C 31.36316 2.30620 -0.54645
C 31.46455 1.27593 -1.64584
C 31.72951 -0.11700 -0.90617
C 31.64282 0.35544 0.56172
C 30.23972 1.14627 -2.55104
C 32.98073 -0.92300 -1.35057
C 34.23841 -0.18345 -0.95474
N 31.61047 0.97197 3.47338
C 31.62076 -0.32756 2.92091
C 31.55626 -1.28034 4.04806
C 31.48918 -0.46234 5.18598
C 31.49000 0.90376 4.79448
C 31.51606 -2.69885 3.93699
C 31.40714 -0.44394 6.62877
O 31.45914 -1.29809 7.51278
C 31.06505 1.03608 7.09664
C 32.00167 1.39728 8.21682
O 33.22350 1.50222 8.01007
O 31.28149 1.77005 9.29048
C 32.16640 2.14840 10.42433
H 29.72671 6.79921 3.46839
H 31.50630 3.84743 -1.74041
H 31.76414 -1.60201 1.36954
H 29.81904 3.89189 7.45753
H 29.89112 6.30496 6.00529
H 27.93694 4.04330 6.08389
H 27.98899 4.92934 4.54554
H 27.61220 5.81547 5.99207
H 32.11142 5.68798 6.62999
H 32.49687 4.15607 7.49593
H 30.49251 5.52890 8.98098
H 31.46574 6.77735 8.28884
H 30.31691 8.40838 2.24243
H 30.08748 8.58978 0.53549
H 31.72459 8.84503 1.18471
H 31.17272 8.60572 -1.41319
H 32.51540 8.30713 -2.66010
H 32.84302 8.45165 -0.90933
H 32.29994 1.49041 -2.30099
H 30.90026 -0.82340 -1.05843
H 30.47045 1.29940 -3.62220
H 29.49264 1.85441 -2.21796
H 29.70902 0.19860 -2.47702
H 33.00465 -1.01547 -2.43917
H 33.01551 -1.92710 -0.94124
H 34.36778 0.53649 -1.76839
H 35.08244 -0.86764 -0.89518
H 34.08544 0.38687 -0.03117
H 32.19015 -3.08917 4.68783
H 31.79806 -3.02765 2.93893
H 30.45842 -2.88092 4.14562
H 30.02215 1.04597 7.41575
H 31.87116 1.53239 11.27935
H 32.06104 3.24425 10.53273
H 33.21909 1.90463 10.34149

