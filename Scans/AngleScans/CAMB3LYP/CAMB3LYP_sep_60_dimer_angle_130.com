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
Mg 31.54051 2.83181 2.39494
C 31.66286 2.13606 -1.12358
C 30.91945 -0.48717 2.87787
C 30.98091 3.48266 5.44956
C 31.03229 6.18046 1.57398
N 31.34102 1.14275 0.98089
C 31.46843 0.98750 -0.32831
C 31.46285 -0.48996 -0.80436
C 30.66244 -1.13680 0.42845
C 31.02000 -0.12133 1.51180
C 29.15454 -1.12353 0.24034
C 32.90990 -1.08452 -0.92504
C 32.99797 -2.41487 -1.67617
H 34.17293 -2.66975 -2.50500
N 31.51164 1.63590 3.97847
C 31.26476 0.30920 3.99182
C 31.34005 -0.09986 5.40291
C 31.48853 1.12269 6.23841
C 31.30794 2.15178 5.23420
C 31.35043 -1.53778 5.83385
C 31.73747 1.33112 7.70155
O 31.59316 2.47949 8.25224
C 32.06430 0.13283 8.48818
N 30.86432 4.53773 3.37216
C 30.75978 4.57115 4.71227
C 30.52074 5.97815 5.20583
C 30.72837 6.92162 3.93134
C 30.97470 5.84199 2.85471
C 29.14512 6.24407 5.81686
C 31.75882 8.07553 4.07031
C 33.15470 7.50869 4.19517
N 31.55662 4.00882 0.58319
C 31.29279 5.39072 0.46140
C 31.30879 5.71157 -0.98055
C 31.56002 4.47994 -1.60381
C 31.67126 3.45967 -0.62087
C 31.06048 6.99197 -1.55068
C 31.74584 3.80906 -2.87034
O 31.84408 4.17485 -4.04093
C 31.69504 2.23976 -2.62319
C 32.85997 1.61813 -3.34364
O 34.02645 1.86505 -2.99015
O 32.40396 0.67182 -4.18443
C 33.52329 0.01871 -4.91407
H 30.49655 -1.46151 3.12427
H 30.88656 3.78162 6.49608
H 30.97966 7.21949 1.24183
H 30.92735 -0.66716 -1.72935
H 31.05385 -2.11858 0.68197
H 28.86338 -0.57257 -0.64160
H 28.75134 -0.65573 1.12895
H 28.76625 -2.14605 0.23544
H 33.24861 -1.40669 0.06804
H 33.57814 -0.37327 -1.37023
H 32.07939 -2.63512 -2.23353
H 33.06739 -3.22290 -0.94386
H 31.06202 -2.20582 5.01291
H 30.55118 -1.65568 6.57354
H 32.29843 -1.83470 6.28750
H 31.25488 -0.58809 8.40124
H 32.29555 0.49596 9.49144
H 32.95344 -0.33848 8.05467
H 31.24354 6.24954 5.96544
H 29.79871 7.43842 3.65119
H 29.19480 6.63077 6.85227
H 28.57242 5.32689 5.78056
H 28.51547 6.93199 5.25517
H 31.57116 8.64269 4.98534
H 31.73370 8.77883 3.24464
H 33.22811 7.26573 5.25944
H 33.89668 8.25088 3.90766
H 33.25016 6.57090 3.63562
H 31.80219 7.13724 -2.32479
H 31.10928 7.77779 -0.79975
H 30.04401 6.84482 -1.92539
H 30.73864 1.87899 -3.00363
H 33.30969 0.11900 -5.98272
H 33.58636 -1.00804 -4.50719
H 34.50120 0.48129 -4.84953

