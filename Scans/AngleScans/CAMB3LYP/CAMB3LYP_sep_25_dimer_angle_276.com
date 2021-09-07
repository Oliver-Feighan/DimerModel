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
Mg 13.14542 1.19033 0.99502
C 12.52266 3.71654 3.46676
C 12.23486 3.48255 -1.36108
C 13.24958 -1.13351 -1.16283
C 12.89509 -1.17085 3.54608
N 12.48474 3.29506 1.15593
C 12.34521 4.16794 2.14230
C 12.07454 5.62446 1.67876
C 11.44798 5.30485 0.23510
C 12.12057 3.95965 -0.03088
C 9.94405 5.08904 0.26615
C 13.39202 6.46823 1.55912
C 13.17866 7.97640 1.41201
H 14.13866 8.87974 2.04039
N 13.27142 1.27972 -0.98354
C 12.87512 2.29526 -1.77909
C 13.16415 1.86431 -3.15565
C 13.61111 0.44531 -3.11887
C 13.36952 0.12888 -1.72546
C 13.08366 2.79247 -4.33271
C 14.15277 -0.47837 -4.16750
O 14.25206 -1.74041 -3.96692
C 14.47603 0.12475 -5.46875
N 12.87244 -0.86667 1.12385
C 13.02554 -1.65333 0.04419
C 13.05258 -3.11255 0.43192
C 13.12868 -3.13133 2.02924
C 13.03853 -1.61430 2.30473
C 11.85745 -3.94591 -0.03022
C 14.29695 -3.93057 2.66880
C 15.61384 -3.25496 2.36102
N 12.96196 1.23934 3.14723
C 12.84655 0.13434 4.01885
C 12.63026 0.67233 5.37760
C 12.61170 2.06360 5.19732
C 12.78274 2.37040 3.82033
C 12.43390 -0.10041 6.55677
C 12.47568 3.34213 5.85729
O 12.39576 3.71113 7.02827
C 12.28546 4.46531 4.74883
C 13.21227 5.60361 5.07716
O 14.44529 5.44790 5.03200
O 12.49738 6.73908 5.17721
C 13.37426 7.89859 5.49097
H 11.75257 4.04857 -2.15851
H 13.38925 -1.97120 -1.85003
H 12.90535 -1.83934 4.40961
H 11.35800 6.16775 2.28307
H 11.76014 6.03967 -0.50243
H 9.56060 5.06948 1.27547
H 9.77313 4.13850 -0.22200
H 9.44271 5.83957 -0.35166
H 13.86993 6.25117 0.59516
H 14.04493 6.27479 2.38792
H 12.15235 8.27360 1.65962
H 13.28726 8.24099 0.35739
H 12.56906 3.72617 -4.07427
H 12.42888 2.31677 -5.07084
H 14.05761 2.97688 -4.79107
H 13.58462 0.59024 -5.88246
H 14.93272 -0.67438 -6.05586
H 15.20560 0.92691 -5.31099
H 13.93198 -3.60261 0.03242
H 12.23133 -3.58546 2.47460
H 12.14594 -4.81898 -0.64551
H 11.18296 -3.30106 -0.57772
H 11.22025 -4.32199 0.76837
H 14.35297 -4.93192 2.23487
H 14.20158 -4.04991 3.74289
H 15.85584 -3.62962 1.36194
H 16.37187 -3.54519 3.08572
H 15.49035 -2.16842 2.28463
H 13.02939 0.36032 7.33370
H 12.71579 -1.13988 6.40245
H 11.35507 0.01697 6.68949
H 11.23867 4.77098 4.76844
H 12.97775 8.36240 6.39946
H 13.38964 8.51811 4.57461
H 14.39497 7.68896 5.78872

