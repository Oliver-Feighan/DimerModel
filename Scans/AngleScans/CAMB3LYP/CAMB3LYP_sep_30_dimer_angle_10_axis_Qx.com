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
Mg 15.77169 1.42553 1.20441
C 16.01937 -0.74929 4.04832
C 14.11280 3.62337 3.21741
C 15.14214 3.16060 -1.37641
C 16.40457 -1.35341 -0.80124
N 15.20123 1.33573 3.33953
C 15.40282 0.48633 4.33564
C 14.98871 1.02281 5.73216
C 13.91049 2.11525 5.26005
C 14.45531 2.40887 3.86376
C 12.50644 1.55197 5.11620
C 16.18630 1.69022 6.49527
C 15.93326 1.95789 7.98056
H 17.04795 1.81582 8.91318
N 15.21954 3.32437 1.03678
C 14.57543 4.06903 1.95976
C 14.38152 5.39167 1.34560
C 14.81680 5.31731 -0.07556
C 15.06226 3.89469 -0.20198
C 13.90538 6.58639 2.11978
C 14.97252 6.35307 -1.14755
O 15.13596 6.03022 -2.37714
C 14.83467 7.75566 -0.72894
N 15.56186 0.93789 -0.80570
C 15.34128 1.89036 -1.72877
C 15.50049 1.34233 -3.12688
C 16.11346 -0.12371 -2.94645
C 16.11954 -0.20742 -1.40426
C 14.22082 1.25979 -3.95876
C 17.43064 -0.43452 -3.70876
C 18.56418 0.38555 -3.13610
N 16.32840 -0.63594 1.53372
C 16.51702 -1.63546 0.55419
C 16.77314 -2.90271 1.26897
C 16.69232 -2.55432 2.62565
C 16.38607 -1.17200 2.74771
C 16.98929 -4.17158 0.66123
C 16.78594 -3.04385 3.98223
O 17.10700 -4.11288 4.49976
C 16.22949 -1.92381 4.96318
C 17.21088 -1.76581 6.09209
O 18.35587 -1.33024 5.87789
O 16.57036 -1.94707 7.26138
C 17.49980 -1.78942 8.41165
H 13.38788 4.28060 3.69848
H 15.04202 3.73586 -2.29974
H 16.70797 -2.24463 -1.35473
H 14.51850 0.29150 6.37858
H 13.95303 3.00624 5.88108
H 12.48777 0.47955 5.24153
H 12.18173 1.82215 4.11990
H 11.82416 2.06008 5.80370
H 16.30959 2.72038 6.13705
H 17.08178 1.11104 6.37990
H 15.05113 1.42201 8.35153
H 15.67759 3.01303 8.10438
H 13.50730 6.29772 3.10040
H 13.03997 6.99069 1.58363
H 14.66605 7.36572 2.20167
H 13.85505 7.90420 -0.28092
H 15.06665 8.34931 -1.61530
H 15.57323 7.96132 0.05391
H 16.19325 1.94458 -3.70175
H 15.42088 -0.90135 -3.30064
H 14.28512 1.81328 -4.91458
H 13.39964 1.63000 -3.35947
H 13.89259 0.25075 -4.20227
H 17.33790 -0.14543 -4.75852
H 17.70511 -1.48389 -3.68675
H 18.45388 1.34802 -3.64459
H 19.52343 -0.07423 -3.36519
H 18.42035 0.56095 -2.06363
H 17.81209 -4.63398 1.19008
H 17.20388 -4.07237 -0.40078
H 16.01909 -4.64211 0.84177
H 15.25081 -2.25288 5.31456
H 17.43390 -2.70742 9.00379
H 17.20270 -0.84691 8.90880
H 18.56135 -1.75210 8.19679

