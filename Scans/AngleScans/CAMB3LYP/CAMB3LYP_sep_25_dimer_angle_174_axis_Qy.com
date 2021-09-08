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
Mg 13.13239 1.18503 0.98972
C 13.99565 -1.50575 -1.22241
C 14.22866 3.32631 -1.42848
C 12.80590 3.55210 3.07726
C 13.18746 -1.12597 3.59778
N 13.98001 0.83499 -1.02305
C 14.20385 -0.22762 -1.78141
C 14.60075 0.09951 -3.24599
C 15.18638 1.57643 -3.01227
C 14.39421 1.95871 -1.76367
C 16.66511 1.58905 -2.66221
C 13.36337 0.11825 -4.21064
C 13.70997 0.12463 -5.70114
H 12.83902 -0.60137 -6.62138
N 13.00152 3.14398 0.69998
C 13.48180 3.84541 -0.34817
C 13.14603 5.25210 -0.07866
C 12.57401 5.34349 1.29216
C 12.79571 3.99032 1.76105
C 13.30136 6.33673 -1.10476
C 11.94471 6.46674 2.05930
O 11.73413 6.38549 3.32100
C 11.66797 7.69968 1.30780
N 13.22088 1.25645 3.06564
C 12.99063 2.40350 3.72835
C 12.83561 2.15583 5.20991
C 12.76894 0.56654 5.37515
C 12.99648 0.14956 3.90554
C 13.94813 2.71330 6.09755
C 11.53813 -0.01199 6.12554
C 10.28500 0.21001 5.30960
N 13.33406 -0.95903 1.16401
C 13.35592 -1.71984 2.35345
C 13.62872 -3.12009 1.96967
C 13.77061 -3.07257 0.57465
C 13.61842 -1.73374 0.12312
C 13.76304 -4.21735 2.86639
C 14.02507 -3.84888 -0.61753
O 14.14564 -5.04826 -0.86386
C 14.30763 -2.84939 -1.82068
C 13.48875 -3.29855 -2.99975
O 12.24646 -3.25768 -2.95884
O 14.30313 -3.49509 -4.05261
C 13.53579 -3.93100 -5.24968
H 14.75430 4.07364 -2.02353
H 12.58710 4.31357 3.82944
H 13.12321 -1.92213 4.34254
H 15.36715 -0.54265 -3.66319
H 14.93632 2.23593 -3.83925
H 17.05209 0.59215 -2.51197
H 16.74689 2.16788 -1.75154
H 17.22748 2.14035 -3.42125
H 12.86143 1.09104 -4.13030
H 12.70138 -0.69832 -3.99708
H 14.76043 -0.13437 -5.88117
H 13.61838 1.14754 -6.07434
H 13.89927 5.99865 -1.96010
H 13.90589 7.12666 -0.64605
H 12.34477 6.76053 -1.41781
H 12.59468 8.08089 0.88548
H 11.13758 8.35301 2.00329
H 11.01427 7.45526 0.46308
H 11.91316 2.58658 5.57957
H 13.62501 0.18014 5.94759
H 13.57844 3.40429 6.87854
H 14.67396 3.20731 5.46536
H 14.55447 1.96391 6.60348
H 11.38971 0.51422 7.07170
H 11.62970 -1.06829 6.35543
H 10.00366 1.23634 5.56371
H 9.50894 -0.49535 5.59983
H 10.50480 0.18468 4.23596
H 13.21648 -5.04356 2.43148
H 13.38814 -3.96928 3.85719
H 14.84892 -4.34409 2.85825
H 15.37769 -2.88193 -2.02919
H 13.97839 -4.87327 -5.58720
H 13.56977 -3.07800 -5.95323
H 12.50244 -4.22307 -5.10360

