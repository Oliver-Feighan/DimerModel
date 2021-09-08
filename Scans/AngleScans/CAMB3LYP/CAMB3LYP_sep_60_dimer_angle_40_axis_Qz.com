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
Mg 31.53925 2.83859 2.40354
C 31.59140 -0.67225 3.14540
C 30.37066 3.33798 5.56892
C 31.05465 5.87947 1.63828
C 31.60366 1.98430 -0.97458
N 31.08155 1.43487 4.05046
C 31.19886 0.12894 4.23785
C 30.95496 -0.33569 5.69874
C 30.04216 0.89014 6.19172
C 30.54856 1.97075 5.23856
C 28.55992 0.67967 5.93095
C 32.28501 -0.43038 6.52579
C 32.16148 -1.16979 7.85986
H 33.28891 -1.97917 8.31405
N 31.29165 4.43073 3.56215
C 30.82828 4.45071 4.82949
C 30.81602 5.86588 5.23111
C 31.15383 6.69399 4.04164
C 31.15955 5.67929 3.00703
C 30.58236 6.30805 6.64650
C 31.41440 8.15899 3.86255
O 31.45499 8.69858 2.70064
C 31.52776 8.95962 5.09041
N 31.14204 3.79247 0.59947
C 31.02668 5.13059 0.53569
C 31.01746 5.60966 -0.89637
C 31.39541 4.33110 -1.77956
C 31.47379 3.26662 -0.66325
C 29.69690 6.19822 -1.39245
C 32.60082 4.47634 -2.74835
C 33.88175 4.62619 -1.95958
N 31.77416 1.01821 1.26376
C 31.74451 0.88183 -0.14144
C 31.83259 -0.56217 -0.44067
C 31.88460 -1.17209 0.82178
C 31.81224 -0.17973 1.83643
C 31.80739 -1.14579 -1.73869
C 31.97354 -2.43050 1.52690
O 32.14656 -3.60231 1.19424
C 31.66023 -2.17198 3.06358
C 32.71557 -2.87031 3.87667
O 33.90204 -2.50154 3.82267
O 32.12033 -3.71038 4.74289
C 33.12561 -4.41833 5.57958
H 29.78898 3.58565 6.45728
H 30.99722 6.92214 1.31729
H 31.72827 1.64337 -2.00459
H 30.40985 -1.26709 5.79418
H 30.26215 1.15701 7.22213
H 28.37579 -0.21071 5.34832
H 28.22801 1.55858 5.39393
H 28.00773 0.67696 6.87510
H 32.55246 0.57007 6.88954
H 33.06769 -0.87114 5.93947
H 31.22660 -1.73894 7.93066
H 32.08637 -0.43031 8.66082
H 30.19821 5.48814 7.26584
H 29.76486 7.03675 6.62303
H 31.46197 6.77794 7.09165
H 30.61133 8.86631 5.66826
H 31.80260 9.96335 4.76046
H 32.33222 8.54295 5.70672
H 31.76507 6.37774 -1.05203
H 30.56800 4.03325 -2.44025
H 29.79618 7.23122 -1.77601
H 28.98072 6.16056 -0.58239
H 29.19740 5.62197 -2.16940
H 32.49755 5.38407 -3.34795
H 32.70351 3.64496 -3.43774
H 33.89975 5.69327 -1.71857
H 34.74016 4.34393 -2.56579
H 33.82798 4.07539 -1.01333
H 32.57313 -1.90997 -1.75138
H 31.97569 -0.40034 -2.51312
H 30.78566 -1.53433 -1.75804
H 30.66243 -2.56370 3.26489
H 32.94579 -5.49076 5.45610
H 33.01233 -4.00262 6.59845
H 34.16568 -4.34294 5.28461

