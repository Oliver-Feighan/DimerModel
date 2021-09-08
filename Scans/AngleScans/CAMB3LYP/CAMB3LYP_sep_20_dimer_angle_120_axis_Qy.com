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
Mg 10.51367 0.94903 0.78861
C 13.12483 -1.42583 0.13980
C 12.83143 3.40450 0.30560
C 8.35526 3.04656 1.79313
C 8.73531 -1.63175 2.31266
N 12.66749 0.87353 0.29100
C 13.53787 -0.08371 0.00729
C 14.90534 0.42637 -0.52130
C 14.87713 1.89837 0.12002
C 13.36548 2.09420 0.21611
C 15.45620 1.95439 1.52390
C 14.95673 0.48350 -2.08850
C 16.35517 0.67541 -2.67950
H 16.67611 0.00912 -3.93868
N 10.42903 2.91770 0.55249
C 11.46452 3.75736 0.34272
C 10.87923 5.10427 0.25658
C 9.43351 5.01013 0.59679
C 9.35262 3.62759 1.02357
C 11.66046 6.30355 -0.19600
C 8.31230 6.00391 0.55814
O 7.18613 5.77083 1.12403
C 8.60291 7.29561 -0.08105
N 8.88900 0.79330 2.07645
C 8.08131 1.84462 2.30069
C 6.83141 1.42498 3.03691
C 6.85543 -0.17405 3.04749
C 8.21981 -0.41116 2.36364
C 6.69818 1.94119 4.46945
C 5.60657 -0.90045 2.47707
C 5.50432 -0.65912 0.98822
N 10.75510 -1.18647 1.01082
C 9.90607 -2.07278 1.70916
C 10.54567 -3.40407 1.67752
C 11.74305 -3.19341 0.97725
C 11.85226 -1.82316 0.61668
C 10.03872 -4.58528 2.28896
C 12.94450 -3.81649 0.47031
O 13.36025 -4.97230 0.39973
C 13.95258 -2.67483 0.01556
C 14.47769 -3.03468 -1.34727
O 13.71580 -3.06803 -2.32962
O 15.82221 -3.06680 -1.30706
C 16.39018 -3.40862 -2.63837
H 13.52376 4.24170 0.39864
H 7.52985 3.70598 2.07123
H 8.19763 -2.50834 2.68028
H 15.76614 -0.12147 -0.15705
H 15.31448 2.63100 -0.55328
H 15.68375 0.97011 1.90542
H 14.70114 2.43169 2.13472
H 16.32579 2.61743 1.54680
H 14.47993 1.41172 -2.42879
H 14.49988 -0.38746 -2.51669
H 17.14386 0.49715 -1.93847
H 16.47566 1.72685 -2.95149
H 12.73758 6.09687 -0.21858
H 11.54710 7.06989 0.57842
H 11.30225 6.70557 -1.14613
H 9.43530 7.77266 0.43081
H 7.65473 7.83666 -0.09092
H 8.93053 7.11083 -1.11020
H 5.94386 1.75965 2.51407
H 6.94210 -0.57350 4.06864
H 5.77036 2.51904 4.64065
H 7.56820 2.54236 4.69793
H 6.73741 1.17506 5.24192
H 4.69531 -0.49216 2.92098
H 5.60529 -1.96903 2.66450
H 5.00992 0.31516 0.92953
H 4.90572 -1.43465 0.51466
H 6.49792 -0.55223 0.53769
H 10.17121 -5.38704 1.57473
H 8.99384 -4.47059 2.56977
H 10.69371 -4.64953 3.16196
H 14.74771 -2.62415 0.76039
H 17.03521 -4.28119 -2.49605
H 16.87012 -2.48191 -3.00528
H 15.70656 -3.77237 -3.39652

