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
Mg 31.52392 2.83882 2.40147
C 32.44343 5.29400 4.85212
C 31.24812 0.60310 4.96274
C 31.15226 0.55203 0.23356
C 32.95731 4.90914 -0.00786
N 31.83100 3.04206 4.58245
C 32.11925 4.03819 5.40645
C 31.97827 3.68860 6.91226
C 32.16877 2.09743 6.80857
C 31.68770 1.88660 5.37442
C 33.62105 1.65831 6.89463
C 30.56507 4.06195 7.48279
C 30.46028 4.03538 9.00925
H 29.60408 5.01970 9.66538
N 30.77431 1.00991 2.57957
C 30.72284 0.25697 3.69838
C 30.09758 -1.01361 3.29982
C 29.93438 -1.01241 1.82061
C 30.64858 0.20099 1.47759
C 29.64702 -2.04545 4.29266
C 29.25959 -1.95495 0.87062
O 29.44549 -1.88260 -0.39549
C 28.44866 -3.02343 1.47247
N 32.17158 2.64066 0.43587
C 31.82682 1.57012 -0.30102
C 32.17234 1.77566 -1.75662
C 32.60217 3.31083 -1.88352
C 32.51239 3.72242 -0.39761
C 33.28812 0.88784 -2.30727
C 31.84647 4.16864 -2.93507
C 30.40400 4.34376 -2.51793
N 32.35764 4.83163 2.36004
C 32.92640 5.49540 1.25100
C 33.45980 6.78302 1.74085
C 33.18263 6.76880 3.11616
C 32.54113 5.54904 3.46282
C 34.14313 7.75336 0.95497
C 33.29639 7.50370 4.35539
O 33.67276 8.63373 4.66379
C 32.93636 6.52731 5.55682
C 31.97743 7.24467 6.46696
O 30.83761 7.54722 6.07213
O 32.48226 7.26120 7.71405
C 31.56183 7.95250 8.65572
H 31.35343 -0.22939 5.65898
H 30.95108 -0.15621 -0.57353
H 33.33132 5.65320 -0.71437
H 32.74320 4.11383 7.55088
H 31.52131 1.57566 7.50852
H 34.29730 2.50022 6.89270
H 33.79152 1.03411 6.02724
H 33.76854 1.01208 7.76466
H 29.85520 3.26330 7.23181
H 30.24623 5.01793 7.11529
H 31.44576 4.00229 9.48955
H 29.98483 3.09804 9.30818
H 30.04857 -1.84319 5.29329
H 30.11049 -2.99401 4.00063
H 28.56231 -2.17151 4.30625
H 29.07315 -3.62434 2.12927
H 27.97577 -3.53879 0.63423
H 27.67852 -2.56518 2.10295
H 31.30988 1.59890 -2.38750
H 33.65745 3.41569 -2.17544
H 32.97632 0.28754 -3.18276
H 33.63746 0.24492 -1.51030
H 34.19694 1.41352 -2.59520
H 31.83044 3.65685 -3.90045
H 32.29074 5.14527 -3.09567
H 29.93068 3.42389 -2.87415
H 29.97105 5.22072 -2.99494
H 30.30948 4.36160 -1.42596
H 33.74973 8.71949 1.24189
H 34.00965 7.56890 -0.10895
H 35.17270 7.57486 1.27652
H 33.86848 6.27293 6.06281
H 32.13606 8.75130 9.13516
H 31.15579 7.16019 9.31231
H 30.73652 8.51137 8.23030

