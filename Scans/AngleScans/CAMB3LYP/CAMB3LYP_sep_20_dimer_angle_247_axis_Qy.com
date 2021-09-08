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
Mg 10.50117 0.94942 0.80136
C 8.52438 -1.84690 -0.27193
C 8.63298 2.88430 -1.29643
C 12.50543 3.33967 1.38194
C 12.88904 -1.33837 1.90126
N 8.82313 0.46591 -0.55608
C 8.11692 -0.61943 -0.83473
C 6.85785 -0.36160 -1.70527
C 7.32868 1.00180 -2.41120
C 8.29734 1.50765 -1.34420
C 8.10848 0.77620 -3.69586
C 5.56711 -0.14438 -0.83994
C 4.25263 -0.19501 -1.62182
H 3.07851 -0.77181 -0.97282
N 10.28090 2.90365 0.53481
C 9.46169 3.51775 -0.34443
C 9.68536 4.96036 -0.16276
C 10.82266 5.14373 0.77938
C 11.27030 3.77285 0.92220
C 8.80777 6.00504 -0.78879
C 11.41793 6.35426 1.43243
O 12.55191 6.30893 2.02813
C 10.67874 7.61515 1.27350
N 12.50805 1.00714 1.34121
C 13.12536 2.17678 1.58370
C 14.47737 1.95840 2.21990
C 14.53747 0.40075 2.57726
C 13.18644 -0.04851 1.97870
C 15.68511 2.32992 1.35983
C 14.85282 0.02868 4.05198
C 13.70879 0.44915 4.94629
N 10.62385 -1.19899 0.99457
C 11.72613 -1.95251 1.45421
C 11.37497 -3.37870 1.29583
C 10.09167 -3.35567 0.72915
C 9.68083 -2.01014 0.52852
C 12.21608 -4.48262 1.61209
C 8.99549 -4.16370 0.24524
O 8.73906 -5.36702 0.24102
C 7.98354 -3.22368 -0.54134
C 6.59209 -3.53584 -0.06286
O 6.25846 -3.29527 1.11078
O 5.82554 -3.86182 -1.11951
C 4.43313 -4.16923 -0.69699
H 8.26113 3.53672 -2.08697
H 13.19267 4.12699 1.70024
H 13.54020 -2.11323 2.31145
H 6.66050 -1.11918 -2.45421
H 6.49783 1.69238 -2.53012
H 8.31984 -0.26993 -3.86044
H 9.02840 1.33493 -3.58356
H 7.58203 1.22911 -4.54088
H 5.53949 0.89665 -0.49320
H 5.53123 -0.84343 -0.02712
H 4.38536 -0.62005 -2.62416
H 3.91930 0.82902 -1.80692
H 8.15717 5.57440 -1.55993
H 9.46519 6.68755 -1.33820
H 8.24183 6.57724 -0.05058
H 10.57452 7.84172 0.21508
H 11.21280 8.34592 1.88394
H 9.66541 7.47867 1.66740
H 14.57215 2.53244 3.13346
H 15.32192 -0.11816 2.00698
H 16.35086 3.07202 1.83956
H 15.32595 2.70004 0.40875
H 16.31332 1.49297 1.05970
H 15.73457 0.57261 4.39955
H 15.04805 -1.02840 4.19847
H 13.91576 1.50759 5.13056
H 13.71692 -0.12171 5.87263
H 12.75149 0.38796 4.41575
H 11.59704 -5.21026 2.11986
H 13.05854 -4.17667 2.22886
H 12.52958 -4.78272 0.60863
H 8.10605 -3.42840 -1.60560
H 4.19945 -5.17075 -1.07119
H 3.81484 -3.33336 -1.07488
H 4.24648 -4.29083 0.36359

