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
Mg 10.50068 0.94890 0.79703
C 9.32431 -1.87166 -1.08435
C 9.59770 2.90169 -1.84983
C 11.93358 3.38754 2.23485
C 12.31667 -1.29046 2.75487
N 9.58786 0.46161 -1.15785
C 9.12782 -0.63698 -1.73733
C 8.34591 -0.38980 -3.05532
C 8.99635 1.02168 -3.45970
C 9.39040 1.51646 -2.06952
C 10.25923 0.88543 -4.29391
C 6.80000 -0.26306 -2.81840
C 5.94732 -0.32999 -4.08734
H 4.64236 -0.98144 -4.01546
N 10.30750 2.90300 0.50933
C 9.90715 3.52353 -0.62024
C 9.95276 4.96413 -0.32568
C 10.56996 5.14931 1.01577
C 10.98910 3.79208 1.30272
C 9.36839 6.00069 -1.24082
C 10.76379 6.35153 1.88943
O 11.53778 6.32399 2.91071
C 10.09411 7.58821 1.46087
N 12.08146 1.06117 2.14271
C 12.47191 2.24205 2.65361
C 13.43521 2.04852 3.80034
C 13.42328 0.47987 4.11140
C 12.48158 0.00528 2.98297
C 14.87151 2.50595 3.54702
C 13.10203 0.05409 5.57023
C 11.66521 0.38818 5.90086
N 10.64753 -1.19955 0.97209
C 11.48962 -1.92908 1.83969
C 11.31802 -3.35937 1.51218
C 10.39801 -3.36129 0.45290
C 10.03808 -2.02529 0.12878
C 12.00431 -4.44247 2.13027
C 9.65759 -4.18915 -0.47182
O 9.49382 -5.40021 -0.61414
C 9.02579 -3.25508 -1.59196
C 7.58237 -3.64320 -1.76057
O 6.76899 -3.46939 -0.83602
O 7.35665 -3.95093 -3.05082
C 5.93589 -4.33182 -3.27032
H 9.56176 3.57409 -2.70743
H 12.37600 4.18610 2.83487
H 12.77332 -2.05684 3.38479
H 8.52747 -1.11986 -3.83495
H 8.25820 1.68365 -3.90498
H 10.57764 -0.14303 -4.37783
H 11.01200 1.47432 -3.78638
H 10.11772 1.35505 -5.27154
H 6.57049 0.75930 -2.49152
H 6.46063 -0.99909 -2.11575
H 6.51652 -0.70361 -4.94720
H 5.66850 0.68744 -4.37213
H 9.13187 5.58012 -2.22600
H 10.15835 6.73277 -1.44055
H 8.51189 6.51610 -0.80102
H 10.43723 7.85822 0.46497
H 10.27711 8.31106 2.25829
H 9.01863 7.39418 1.38136
H 13.10110 2.58388 4.68061
H 14.40292 0.01866 3.91794
H 15.22841 3.25122 4.28273
H 14.93078 2.90434 2.54287
H 15.61267 1.70888 3.52345
H 13.72116 0.61605 6.27381
H 13.27430 -1.00003 5.76040
H 11.71574 1.44452 6.18135
H 11.31032 -0.22330 6.72801
H 11.02904 0.31340 5.01128
H 11.26921 -5.21620 2.30764
H 12.48663 -4.13180 3.05472
H 12.73046 -4.68420 1.34954
H 9.59996 -3.40638 -2.50685
H 5.93880 -5.32379 -3.73255
H 5.49185 -3.50471 -3.85553
H 5.32327 -4.50861 -2.39403

