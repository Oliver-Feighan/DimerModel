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
Mg 7.87449 0.71023 0.60493
C 8.93116 3.09007 3.07452
C 7.44254 -1.51683 3.15212
C 7.37338 -1.53748 -1.57776
C 9.45034 2.69861 -1.78437
N 8.18002 0.88254 2.78878
C 8.52516 1.85432 3.61997
C 8.35301 1.50656 5.12296
C 8.44372 -0.09289 5.01162
C 7.95937 -0.26566 3.57349
C 9.86494 -0.62295 5.10389
C 6.96254 1.96515 5.68713
C 6.84670 1.93739 7.21276
H 6.04999 2.97025 7.86921
N 7.01024 -1.06878 0.76846
C 6.90449 -1.82272 1.88276
C 6.20308 -3.04938 1.47349
C 6.04958 -3.03032 -0.00664
C 6.84081 -1.86251 -0.33869
C 5.68230 -4.05591 2.45792
C 5.32283 -3.92364 -0.96581
O 5.52087 -3.85663 -2.23038
C 4.44255 -4.94206 -0.37471
N 8.52074 0.48181 -1.35786
C 8.11398 -0.56111 -2.10268
C 8.48088 -0.37025 -3.55504
C 9.00719 1.13546 -3.67091
C 8.93414 1.54425 -2.18331
C 9.54205 -1.32367 -4.10394
C 8.31355 2.04451 -4.72218
C 6.88236 2.30790 -4.31267
N 8.83211 2.64682 0.57946
C 9.44845 3.27918 -0.52250
C 10.05868 4.52815 -0.02238
C 9.77253 4.52433 1.35116
C 9.05342 3.34559 1.68725
C 10.80661 5.45760 -0.79881
C 9.92449 5.24424 2.59507
O 10.36923 6.34675 2.91194
C 9.49623 4.28628 3.78893
C 8.57861 5.05788 4.69730
O 7.46257 5.43356 4.29738
O 9.07563 5.03620 5.94745
C 8.19460 5.77920 6.88742
H 7.49091 -2.35786 3.84437
H 7.13314 -2.22751 -2.38993
H 9.87483 3.42129 -2.48453
H 9.13915 1.87953 5.76845
H 7.76035 -0.57647 5.70481
H 10.59280 0.17474 5.11064
H 10.00128 -1.25217 4.23410
H 9.96604 -1.28164 5.97121
H 6.20544 1.21405 5.42752
H 6.70676 2.94117 5.32301
H 7.82510 1.83989 7.69873
H 6.31137 1.03030 7.50368
H 6.08948 -3.88446 3.46202
H 6.08702 -5.03025 2.16342
H 4.59174 -4.11355 2.46436
H 5.02388 -5.58444 0.28246
H 3.94346 -5.42234 -1.21858
H 3.69880 -4.43951 0.25370
H 7.61299 -0.48915 -4.19201
H 10.06879 1.17522 -3.95596
H 9.19863 -1.89866 -4.98457
H 9.84525 -1.99138 -3.30847
H 10.48393 -0.85475 -4.38355
H 8.27143 1.53971 -5.69045
H 8.81936 2.99207 -4.87473
H 6.35437 1.42148 -4.67677
H 6.50842 3.21280 -4.78740
H 6.78229 2.32605 -3.22120
H 10.47294 6.44509 -0.50891
H 10.66848 5.28737 -1.86451
H 11.82087 5.21302 -0.47213
H 10.40730 3.97117 4.29906
H 8.81490 6.53781 7.37468
H 7.73541 5.01065 7.53721
H 7.40876 6.39107 6.46020

