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
Mg 31.54058 2.83625 2.38888
C 31.28557 5.43782 -0.06996
C 31.05076 0.60198 -0.14147
C 31.31046 0.55447 4.58158
C 30.87010 5.25068 4.81102
N 31.22223 3.10699 0.21687
C 31.20640 4.13871 -0.61357
C 31.18052 3.75380 -2.11714
C 30.53693 2.28899 -1.97932
C 30.99039 1.95536 -0.55956
C 29.01753 2.29422 -2.00662
C 32.61791 3.69306 -2.74342
C 32.65117 3.62745 -4.27193
H 33.73185 4.31355 -4.97469
N 31.71474 0.86689 2.21615
C 31.50030 0.12498 1.10942
C 31.74135 -1.27054 1.50762
C 31.95489 -1.30802 2.97994
C 31.64083 0.06351 3.32685
C 31.83148 -2.39527 0.51753
C 32.36059 -2.40100 3.92174
O 32.25315 -2.26948 5.19217
C 32.79945 -3.66394 3.31029
N 30.93987 2.84882 4.37919
C 30.98455 1.72885 5.12189
C 30.77091 2.03567 6.58497
C 30.81434 3.63016 6.70036
C 30.96304 3.98812 5.20532
C 29.46637 1.51726 7.18981
C 31.82808 4.23537 7.70968
C 33.24138 3.97749 7.23908
N 31.32872 4.98610 2.42517
C 31.02245 5.79430 3.54190
C 30.87062 7.17863 3.04893
C 31.07804 7.07525 1.66518
C 31.32064 5.71892 1.31730
C 30.53208 8.31149 3.84149
C 31.13642 7.80350 0.41821
O 31.09553 8.99227 0.10375
C 31.14820 6.75702 -0.77809
C 32.23909 7.15620 -1.73368
O 33.43194 7.11388 -1.38450
O 31.71331 7.31218 -2.96234
C 32.75856 7.69753 -3.94761
H 30.67960 -0.16741 -0.81909
H 31.32668 -0.17660 5.39322
H 30.75782 6.07606 5.51728
H 30.54936 4.38054 -2.73570
H 30.97577 1.59605 -2.69249
H 28.61783 3.29734 -1.99545
H 28.70558 1.75279 -1.12311
H 28.65384 1.71418 -2.85962
H 33.07219 2.72302 -2.50380
H 33.21650 4.51590 -2.40393
H 31.68133 3.88163 -4.71661
H 32.81960 2.59009 -4.57115
H 31.46826 -2.09044 -0.47169
H 31.12240 -3.16453 0.84204
H 32.83045 -2.83376 0.46796
H 32.00152 -4.05953 2.68636
H 33.13292 -4.29004 4.14016
H 33.64497 -3.45543 2.64528
H 31.56765 1.61786 7.18808
H 29.84806 4.04148 7.02738
H 29.62252 0.85739 8.06403
H 28.91363 0.99999 6.41686
H 28.76300 2.28796 7.50041
H 31.73102 3.74737 8.68261
H 31.69559 5.30031 7.86891
H 33.43824 2.96153 7.59427
H 33.93011 4.69202 7.68533
H 33.29464 3.96006 6.14431
H 31.17955 9.11811 3.52424
H 30.64682 8.10263 4.90308
H 29.48370 8.44054 3.55949
H 30.16360 6.78383 -1.24651
H 32.42509 8.62654 -4.42029
H 32.88927 6.81694 -4.60430
H 33.72710 7.99264 -3.56123

