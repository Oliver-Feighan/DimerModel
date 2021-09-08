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
Mg 31.52287 2.84196 2.39825
C 29.54609 0.04564 1.32495
C 29.65468 4.77683 0.30046
C 33.52714 5.23220 2.97883
C 33.91074 0.55417 3.49814
N 29.84484 2.35844 1.04080
C 29.13862 1.27310 0.76215
C 27.87956 1.53093 -0.10839
C 28.35038 2.89433 -0.81432
C 29.31904 3.40018 0.25268
C 29.13019 2.66873 -2.09898
C 26.58881 1.74815 0.75695
C 25.27434 1.69752 -0.02494
H 24.10021 1.12072 0.62406
N 31.30260 4.79619 2.13169
C 30.48340 5.41029 1.25246
C 30.70707 6.85290 1.43412
C 31.84437 7.03627 2.37627
C 32.29201 5.66539 2.51909
C 29.82948 7.89758 0.80809
C 32.43963 8.24680 3.02932
O 33.57362 8.20146 3.62501
C 31.70044 9.50768 2.87039
N 33.52975 2.89968 2.93809
C 34.14706 4.06932 3.18058
C 35.49908 3.85093 3.81679
C 35.55917 2.29328 4.17414
C 34.20814 1.84402 3.57558
C 36.70682 4.22246 2.95671
C 35.87452 1.92121 5.64886
C 34.73050 2.34168 6.54318
N 31.64555 0.69355 2.59146
C 32.74783 -0.05998 3.05109
C 32.39667 -1.48617 2.89271
C 31.11338 -1.46314 2.32603
C 30.70253 -0.11761 2.12540
C 33.23778 -2.59009 3.20897
C 30.01720 -2.27117 1.84212
O 29.76077 -3.47448 1.83791
C 29.00525 -1.33115 1.05555
C 27.61380 -1.64331 1.53403
O 27.28016 -1.40273 2.70766
O 26.84725 -1.96928 0.47738
C 25.45484 -2.27669 0.89989
H 29.28283 5.42925 -0.49009
H 34.21437 6.01952 3.29713
H 34.56190 -0.22070 3.90833
H 27.68220 0.77336 -0.85732
H 27.51953 3.58492 -0.93323
H 29.34154 1.62260 -2.26355
H 30.05011 3.22746 -1.98667
H 28.60374 3.12164 -2.94399
H 26.56119 2.78919 1.10369
H 26.55293 1.04910 1.56976
H 25.40706 1.27249 -1.02727
H 24.94101 2.72156 -0.21004
H 29.17887 7.46694 0.03695
H 30.48690 8.58008 0.25868
H 29.26353 8.46978 1.54630
H 31.59622 9.73426 1.81196
H 32.23450 10.23846 3.48082
H 30.68712 9.37120 3.26428
H 35.59385 4.42498 4.73034
H 36.34362 1.77438 3.60387
H 37.37256 4.96455 3.43645
H 36.34766 4.59257 2.00563
H 37.33502 3.38551 2.65659
H 36.75627 2.46515 5.99643
H 36.06976 0.86413 5.79535
H 34.93746 3.40012 6.72744
H 34.73863 1.77083 7.46951
H 33.77320 2.28049 6.01263
H 32.61875 -3.31773 3.71674
H 34.08025 -2.28413 3.82575
H 33.55128 -2.89018 2.20551
H 29.12775 -1.53587 -0.00871
H 25.22116 -3.27822 0.52570
H 24.83655 -1.44082 0.52201
H 25.26819 -2.39830 1.96047

