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
Mg 31.52583 2.83418 2.40209
C 30.92029 5.52269 4.70085
C 32.72646 1.04331 5.04517
C 32.48559 0.62221 0.33969
C 31.36271 5.18082 -0.16937
N 31.74301 3.32756 4.54724
C 31.42193 4.35696 5.31638
C 31.57600 4.09519 6.83865
C 32.66360 2.91561 6.77094
C 32.34112 2.36518 5.38324
C 34.09750 3.41880 6.76382
C 30.23928 3.60678 7.49945
C 30.23293 3.62835 9.02956
H 28.98877 3.96556 9.71574
N 32.00681 0.93067 2.69169
C 32.45505 0.36937 3.83419
C 32.68352 -1.05006 3.52244
C 32.49060 -1.24604 2.05988
C 32.33721 0.12865 1.62769
C 32.96912 -2.07912 4.57738
C 32.46312 -2.46750 1.19185
O 32.51857 -2.38675 -0.08613
C 32.46362 -3.76424 1.88464
N 32.08461 2.91954 0.40138
C 32.40754 1.80291 -0.27446
C 32.50546 2.07181 -1.75719
C 31.94183 3.55420 -1.96313
C 31.68758 3.93543 -0.48820
C 33.90642 1.97194 -2.36029
C 30.78340 3.73137 -2.98261
C 29.53367 3.05530 -2.46638
N 31.02166 4.92863 2.23743
C 31.04366 5.72201 1.06946
C 30.73478 7.10756 1.47844
C 30.57591 7.02766 2.87016
C 30.79183 5.69073 3.30085
C 30.68168 8.23767 0.61463
C 30.28511 7.77231 4.07416
O 29.93505 8.92618 4.31870
C 30.61959 6.85538 5.32865
C 29.46037 6.93427 6.28390
O 28.34639 6.48296 5.96456
O 29.90891 7.32901 7.48945
C 28.79741 7.41161 8.47411
H 33.33076 0.48084 5.75744
H 32.70785 -0.12262 -0.42797
H 31.19676 5.95232 -0.92429
H 31.96846 4.93070 7.40581
H 32.47763 2.16270 7.53246
H 34.14652 4.49476 6.68600
H 34.56754 2.95564 5.90615
H 34.63314 3.04329 7.64041
H 30.12728 2.52872 7.32628
H 29.40338 4.16619 7.12670
H 31.06703 4.21246 9.43703
H 30.41431 2.61347 9.39168
H 33.21475 -1.61190 5.53903
H 33.89031 -2.59287 4.28184
H 32.16904 -2.81607 4.67448
H 33.34857 -3.83803 2.51225
H 32.35171 -4.51488 1.09985
H 31.59802 -3.80290 2.55534
H 31.88815 1.38017 -2.31725
H 32.71928 4.23774 -2.33485
H 33.97301 1.24466 -3.19139
H 34.59992 1.71246 -1.57146
H 34.31778 2.90948 -2.73042
H 31.03270 3.24283 -3.92778
H 30.55936 4.76915 -3.20577
H 29.67966 2.01093 -2.75803
H 28.64783 3.47653 -2.93745
H 29.49162 3.08901 -1.37151
H 29.80644 8.80634 0.89938
H 30.63920 7.93770 -0.43035
H 31.63059 8.71941 0.86516
H 31.54218 7.23136 5.77255
H 28.80926 8.42588 8.88497
H 28.96392 6.57896 9.18314
H 27.78474 7.34954 8.09320

