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
Mg 15.77417 1.41872 1.19063
C 15.41628 4.54476 -0.53529
C 15.25765 -0.12462 -1.80703
C 15.63540 -1.34197 2.74844
C 15.12533 3.13814 4.15175
N 15.39762 2.21518 -0.83819
C 15.34454 3.42071 -1.38429
C 15.28750 3.42217 -2.93552
C 14.67120 1.95557 -3.15486
C 15.16513 1.28862 -1.87272
C 13.15171 1.93432 -3.14951
C 16.70967 3.55076 -3.58582
C 16.70598 3.86888 -5.08267
H 17.75767 4.73185 -5.61343
N 15.97568 -0.44119 0.52828
C 15.74577 -0.88838 -0.72410
C 16.01908 -2.33352 -0.69179
C 16.26972 -2.73206 0.71996
C 15.94235 -1.49743 1.40450
C 16.10267 -3.17378 -1.93292
C 16.71626 -4.01615 1.35077
O 16.63833 -4.20776 2.61570
C 17.16006 -5.07699 0.43472
N 15.22300 0.92179 3.13269
C 15.30415 -0.34668 3.57123
C 15.12203 -0.41888 5.06857
C 15.14265 1.09710 5.57746
C 15.24835 1.81949 4.21644
C 13.84145 -1.09997 5.55087
C 16.17128 1.45363 6.68537
C 17.57642 1.35203 6.13705
N 15.52870 3.48655 1.76669
C 15.23730 3.98408 3.05574
C 15.05100 5.44393 2.92711
C 15.22560 5.69320 1.55748
C 15.48128 4.47198 0.87723
C 14.71407 6.33589 3.98400
C 15.24123 6.71033 0.53088
O 15.17341 7.93878 0.52403
C 15.24012 5.99550 -0.88877
C 16.30035 6.64390 -1.73613
O 17.50204 6.54191 -1.43256
O 15.74176 7.08966 -2.87615
C 16.75585 7.73107 -3.75479
H 14.88220 -0.70876 -2.64768
H 15.68356 -2.25171 3.35138
H 15.01736 3.75883 5.04381
H 14.63115 4.16938 -3.36524
H 15.10328 1.47195 -4.02709
H 12.73631 2.89405 -2.88027
H 12.87057 1.18313 -2.42307
H 12.77631 1.57737 -4.11284
H 17.18532 2.56176 -3.60511
H 17.30319 4.27588 -3.06375
H 15.72141 4.20466 -5.43027
H 16.88360 2.94275 -5.63470
H 15.71012 -2.64004 -2.80723
H 15.41436 -4.01491 -1.79656
H 17.10703 -3.56420 -2.11040
H 16.35335 -5.32191 -0.25208
H 17.52409 -5.88272 1.07516
H 17.98532 -4.69096 -0.17401
H 15.94014 -0.95635 5.53217
H 14.17829 1.39279 6.01612
H 14.02988 -1.95332 6.22938
H 13.27806 -1.42022 4.68451
H 13.13370 -0.44650 6.05808
H 16.10630 0.73657 7.50741
H 16.02566 2.44214 7.10803
H 17.79838 0.28414 6.22335
H 18.26443 1.94764 6.73364
H 17.60273 1.60913 5.07175
H 15.34038 7.21003 3.86525
H 14.85851 5.87161 4.95732
H 13.65707 6.50826 3.76423
H 14.24388 6.11673 -1.31582
H 16.39582 8.74107 -3.97380
H 16.88436 7.04498 -4.61300
H 17.72882 7.94162 -3.32647

