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
Mg 13.13959 1.19142 1.00419
C 14.73815 -0.22395 3.88870
C 10.58800 2.13015 3.06399
C 11.50306 2.04592 -1.57632
C 15.20648 -0.81752 -0.95485
N 12.81981 0.93052 3.17719
C 13.52404 0.43138 4.18181
C 12.93087 0.71708 5.58755
C 11.39767 0.93413 5.16202
C 11.60685 1.39505 3.72102
C 10.59133 -0.35375 5.13925
C 13.53463 2.01031 6.23951
C 13.23352 2.17960 7.73029
H 14.25462 2.78248 8.58254
N 11.56775 2.38865 0.81839
C 10.64692 2.67507 1.76248
C 9.68551 3.58680 1.12326
C 10.02227 3.68479 -0.32301
C 11.05389 2.67192 -0.42267
C 8.62870 4.32478 1.89266
C 9.49328 4.53869 -1.43532
O 9.76512 4.28982 -2.66292
C 8.57224 5.61850 -1.05168
N 13.17554 0.53722 -0.96890
C 12.39822 1.11340 -1.90255
C 12.79252 0.66878 -3.29080
C 14.15876 -0.14261 -3.11021
C 14.27618 -0.10105 -1.57077
C 11.77489 -0.20525 -4.02352
C 15.37330 0.32685 -3.95708
C 15.82759 1.69268 -3.49503
N 14.81763 -0.12333 1.35588
C 15.51899 -0.88645 0.39681
C 16.50205 -1.71009 1.13005
C 16.28700 -1.38346 2.47747
C 15.22992 -0.43907 2.57850
C 17.39967 -2.64901 0.54796
C 16.70670 -1.63080 3.83812
O 17.61719 -2.26988 4.36378
C 15.63762 -0.98604 4.82186
C 16.38232 -0.20570 5.87020
O 17.04033 0.80267 5.55905
O 16.02041 -0.64761 7.08844
C 16.72428 0.10343 8.16190
H 9.63544 2.26838 3.57623
H 11.04531 2.38828 -2.50730
H 15.95403 -1.39676 -1.50105
H 13.00924 -0.10466 6.28924
H 10.93207 1.72068 5.75016
H 11.21372 -1.22165 5.29901
H 10.12930 -0.39430 4.16153
H 9.76953 -0.29674 5.85880
H 13.01202 2.88952 5.84124
H 14.59373 2.06013 6.07684
H 12.85310 1.25515 8.18154
H 12.41022 2.88971 7.83956
H 8.51795 3.92546 2.90837
H 7.67029 4.10685 1.40913
H 8.78614 5.40547 1.89414
H 7.71279 5.19454 -0.53780
H 8.37304 6.17301 -1.97076
H 9.07879 6.27120 -0.33195
H 12.97270 1.52177 -3.93348
H 14.04410 -1.20063 -3.38818
H 11.46127 0.21376 -4.99828
H 10.91875 -0.34703 -3.37741
H 12.09510 -1.22879 -4.21069
H 15.08506 0.43394 -5.00570
H 16.21433 -0.35760 -3.92256
H 15.15030 2.36999 -4.02383
H 16.86308 1.86849 -3.77947
H 15.65207 1.82330 -2.42089
H 18.35767 -2.50346 1.02915
H 17.47078 -2.51577 -0.52952
H 16.90189 -3.58539 0.81381
H 15.05659 -1.80153 5.25422
H 17.23663 -0.63559 8.78563
H 15.94927 0.72390 8.64987
H 17.54977 0.74141 7.86865

