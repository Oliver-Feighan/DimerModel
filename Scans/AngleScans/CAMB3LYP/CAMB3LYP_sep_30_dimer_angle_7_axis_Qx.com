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
Mg 15.77202 1.42515 1.20447
C 15.86457 -0.77553 4.03776
C 14.24201 3.71263 3.21968
C 15.26906 3.20962 -1.37043
C 16.24149 -1.37778 -0.81269
N 15.18364 1.36045 3.33564
C 15.32515 0.49494 4.32822
C 14.93683 1.04923 5.72519
C 13.93243 2.20976 5.25272
C 14.50340 2.47568 3.86134
C 12.49667 1.73669 5.09741
C 16.16920 1.63604 6.49909
C 15.92417 1.91147 7.98430
H 17.02183 1.69476 8.92273
N 15.34144 3.35580 1.04405
C 14.73965 4.13475 1.96728
C 14.63316 5.47011 1.35930
C 15.07183 5.37581 -0.05965
C 15.22813 3.94123 -0.19245
C 14.22824 6.68843 2.13721
C 15.29911 6.40522 -1.12495
O 15.44965 6.07905 -2.35531
C 15.24711 7.81153 -0.69943
N 15.54458 0.96200 -0.80951
C 15.39013 1.93120 -1.72860
C 15.52335 1.38142 -3.12874
C 16.04176 -0.12119 -2.95277
C 16.03288 -0.21304 -1.41106
C 14.24628 1.38383 -3.96866
C 17.34155 -0.51034 -3.70894
C 18.52081 0.23383 -3.12502
N 16.19592 -0.66894 1.52570
C 16.32745 -1.67330 0.54181
C 16.49887 -2.95780 1.25109
C 16.43160 -2.61200 2.60918
C 16.21213 -1.21379 2.73704
C 16.63862 -4.23463 0.63765
C 16.48573 -3.11341 3.96356
O 16.73567 -4.20316 4.47708
C 15.99466 -1.96563 4.94735
C 16.97693 -1.87549 6.08293
O 18.14837 -1.51173 5.87796
O 16.31895 -2.02209 7.24738
C 17.24921 -1.92915 8.40401
H 13.55685 4.41170 3.70004
H 15.21111 3.79477 -2.29114
H 16.49171 -2.28346 -1.36928
H 14.41751 0.34565 6.36476
H 14.02701 3.09311 5.87891
H 12.40981 0.66694 5.21671
H 12.19586 2.03188 4.10070
H 11.84339 2.28318 5.78363
H 16.35928 2.65822 6.14730
H 17.02719 1.00227 6.38586
H 15.00777 1.43026 8.34705
H 15.73458 2.97996 8.11242
H 13.80664 6.42035 3.11383
H 13.39335 7.14913 1.59816
H 15.03588 7.41793 2.22793
H 14.27597 8.01912 -0.25643
H 15.52153 8.39396 -1.58111
H 15.99220 7.96629 0.08892
H 16.25621 1.94183 -3.69615
H 15.30389 -0.85189 -3.31536
H 14.35126 1.93708 -4.92101
H 13.44626 1.80189 -3.37223
H 13.85677 0.39870 -4.21969
H 17.27378 -0.21060 -4.75761
H 17.54935 -1.57500 -3.69109
H 18.47444 1.20393 -3.62884
H 19.45067 -0.28421 -3.35094
H 18.38156 0.41242 -2.05247
H 17.42737 -4.75059 1.16883
H 16.86569 -4.14366 -0.42250
H 15.63964 -4.64410 0.80981
H 14.99503 -2.23427 5.29108
H 17.12199 -2.84421 8.99067
H 17.00885 -0.97238 8.90457
H 18.31234 -1.95758 8.19568

