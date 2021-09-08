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
Mg 8.94056 0.81203 0.68227
C 8.49584 0.45641 4.22555
C 7.74686 3.99123 1.00225
C 8.89455 1.00731 -2.48430
C 8.99993 -2.66670 0.48073
N 8.24465 1.96140 2.43923
C 8.17909 1.74377 3.74410
C 7.82729 2.99791 4.58846
C 7.07199 3.84739 3.45410
C 7.75000 3.26101 2.21752
C 5.58424 3.54943 3.36960
C 9.10477 3.74473 5.11031
C 8.83937 4.79370 6.19242
H 9.83226 4.97088 7.24848
N 8.87155 2.40035 -0.50561
C 8.37388 3.61874 -0.20706
C 8.53891 4.42978 -1.42326
C 9.02372 3.54340 -2.51597
C 8.92677 2.25028 -1.86904
C 8.32067 5.91479 -1.44682
C 9.48463 3.81343 -3.91627
O 9.63539 2.86888 -4.76946
C 9.66441 5.22482 -4.28660
N 8.73443 -0.61501 -0.81554
C 8.80220 -0.26985 -2.11330
C 8.90432 -1.48970 -2.99757
C 9.13558 -2.72191 -2.00471
C 9.03356 -1.98076 -0.65352
C 7.69183 -1.77176 -3.88453
C 10.38060 -3.61289 -2.26687
C 11.64356 -2.82810 -1.99408
N 8.96569 -0.82728 2.08943
C 8.96393 -2.20754 1.79117
C 8.86678 -2.93044 3.07585
C 8.79499 -1.91224 4.03845
C 8.82359 -0.64427 3.39740
C 8.80578 -4.34390 3.23248
C 8.69080 -1.62178 5.45035
O 8.71546 -2.29501 6.47986
C 8.36502 -0.07596 5.62558
C 9.28971 0.47757 6.67490
O 10.51611 0.52627 6.47448
O 8.55891 1.05472 7.64608
C 9.43195 1.62831 8.70467
H 7.17502 4.91843 0.95557
H 8.98838 1.01738 -3.57268
H 9.11119 -3.75273 0.50967
H 7.15905 2.81232 5.42080
H 7.29205 4.90809 3.54273
H 5.30120 2.72652 4.00902
H 5.39126 3.30669 2.33283
H 5.00612 4.45377 3.58004
H 9.49228 4.39152 4.31266
H 9.83982 3.04505 5.45783
H 7.83464 4.69797 6.62171
H 8.83823 5.78153 5.72530
H 7.80971 6.26071 -0.53982
H 7.60979 6.12073 -2.25435
H 9.24044 6.47276 -1.63483
H 8.72535 5.75614 -4.15119
H 10.08237 5.21138 -5.29506
H 10.38495 5.67843 -3.59691
H 9.75311 -1.41287 -3.66596
H 8.29769 -3.43377 -2.03670
H 7.94161 -1.82575 -4.96099
H 6.95099 -1.00421 -3.70410
H 7.14506 -2.68180 -3.64354
H 10.42000 -3.90874 -3.31810
H 10.39305 -4.52245 -1.67566
H 11.79700 -2.27573 -2.92597
H 12.47559 -3.49952 -1.79165
H 11.48531 -2.09236 -1.19702
H 9.46182 -4.59455 4.05546
H 9.09803 -4.85809 2.31926
H 7.74227 -4.47007 3.45222
H 7.31756 0.00698 5.91803
H 9.11370 1.18790 9.65464
H 9.34235 2.72624 8.60337
H 10.48181 1.35999 8.68835

