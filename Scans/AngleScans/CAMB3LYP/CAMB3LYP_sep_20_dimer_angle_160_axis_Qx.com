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
Mg 10.50236 0.94528 0.80470
C 11.13438 3.52423 3.21896
C 10.53915 -1.27803 3.39124
C 10.39591 -1.39466 -1.33558
C 11.63642 3.15148 -1.64321
N 10.80661 1.21030 2.97943
C 10.97693 2.24407 3.78978
C 10.89836 1.89674 5.30066
C 11.28603 0.34098 5.21023
C 10.81874 0.05511 3.78454
C 12.78286 0.08858 5.28190
C 9.45614 2.09628 5.88563
C 9.37314 2.07419 7.41350
H 8.40770 2.95071 8.07085
N 9.99053 -0.96110 1.01048
C 10.04696 -1.70173 2.13722
C 9.58170 -3.04519 1.75931
C 9.40259 -3.08135 0.28238
C 9.95470 -1.79191 -0.08176
C 9.27582 -4.11402 2.76800
C 8.84063 -4.11189 -0.64949
O 9.00135 -4.03124 -1.91854
C 8.17734 -5.26676 -0.02673
N 11.14706 0.80755 -1.16703
C 10.93103 -0.30613 -1.88882
C 11.23117 -0.07547 -3.35070
C 11.46326 1.49995 -3.49827
C 11.33964 1.91394 -2.01547
C 12.44321 -0.82248 -3.90689
C 10.59371 2.24410 -4.54842
C 9.14558 2.24157 -4.11448
N 11.07860 3.02634 0.73287
C 11.54664 3.74355 -0.38990
C 11.91968 5.09343 0.08032
C 11.66238 5.06022 1.45915
C 11.18323 3.77364 1.82606
C 12.46660 6.13282 -0.72384
C 11.69721 5.81767 2.68950
O 11.93218 6.98945 2.98186
C 11.47658 4.81759 3.90504
C 10.44569 5.41929 4.82033
O 9.27242 5.57189 4.43726
O 10.95880 5.51325 6.06060
C 9.96975 6.09429 7.00717
H 10.75623 -2.08272 4.09425
H 10.27600 -2.13171 -2.13299
H 11.90583 3.92849 -2.36178
H 11.61114 2.42185 5.92514
H 10.71736 -0.24989 5.92357
H 13.34793 1.00863 5.26296
H 13.02038 -0.51914 4.41859
H 13.02039 -0.52406 6.15630
H 8.84937 1.21206 5.65170
H 9.01552 3.00044 5.51284
H 10.36044 2.17050 7.88137
H 9.02268 1.08807 7.72776
H 9.66030 -3.85157 3.76133
H 9.85137 -5.00020 2.47940
H 8.21577 -4.37505 2.79680
H 8.87991 -5.77702 0.62786
H 7.76330 -5.84692 -0.85369
H 7.36306 -4.90161 0.60906
H 10.39053 -0.36624 -3.96866
H 12.49357 1.73307 -3.80478
H 12.19921 -1.46711 -4.77233
H 12.87971 -1.40735 -3.10819
H 13.27542 -0.19023 -4.21167
H 10.63095 1.72336 -5.50842
H 10.90990 3.26687 -4.72435
H 8.78748 1.26557 -4.45548
H 8.60042 3.05177 -4.59451
H 9.06217 2.25987 -3.02162
H 11.95829 7.04516 -0.44142
H 12.34508 5.92095 -1.78408
H 13.51407 6.08860 -0.41382
H 12.43904 4.68799 4.40148
H 10.44458 6.96425 7.47128
H 9.67405 5.26483 7.67673
H 9.07592 6.54029 6.58693

