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
Mg 8.93281 0.79695 0.68615
C 6.74737 2.53056 2.94392
C 11.14883 0.53836 3.26628
C 11.05837 -0.32896 -1.38308
C 7.09120 2.17300 -1.93313
N 8.84500 1.48451 2.78659
C 7.93864 2.07806 3.54852
C 8.29644 2.12152 5.05838
C 9.89271 2.00100 4.93078
C 9.97625 1.26411 3.59561
C 10.58745 3.34404 4.77938
C 7.69054 0.91068 5.85133
C 7.73771 1.05116 7.37442
H 6.63317 0.50527 8.15829
N 10.60150 -0.23302 0.99144
C 11.36192 -0.23849 2.10635
C 12.48880 -1.13988 1.82034
C 12.41971 -1.53727 0.38785
C 11.34567 -0.67994 -0.07204
C 13.44858 -1.60341 2.87740
C 13.20071 -2.51643 -0.43525
O 13.12987 -2.52782 -1.71497
C 14.12163 -3.39944 0.29529
N 9.19261 1.07127 -1.35825
C 10.16428 0.42351 -2.02472
C 9.98616 0.56046 -3.51801
C 8.55029 1.23524 -3.71932
C 8.16821 1.46318 -2.24031
C 11.04514 1.39271 -4.24072
C 7.54362 0.48763 -4.63602
C 7.12340 -0.81236 -3.98886
N 7.12162 1.96033 0.50130
C 6.54198 2.45299 -0.68840
C 5.38424 3.28313 -0.29733
C 5.38439 3.23419 1.10485
C 6.47782 2.44682 1.55644
C 4.53212 3.99480 -1.18814
C 4.71445 3.67828 2.30596
O 3.67876 4.29817 2.54433
C 5.64139 3.34636 3.55375
C 4.78757 2.69350 4.60603
O 4.27516 1.57883 4.40216
O 4.89449 3.38794 5.75350
C 4.07403 2.77280 6.83066
H 12.00458 0.60220 3.93895
H 11.69769 -0.78311 -2.14367
H 6.40816 2.55643 -2.69412
H 8.03224 3.04395 5.56146
H 10.30789 1.39158 5.72928
H 9.88082 4.15238 4.66356
H 11.20924 3.25589 3.89811
H 11.27215 3.51068 5.61589
H 8.34198 0.03666 5.72300
H 6.68359 0.71529 5.53748
H 7.95967 2.07880 7.68690
H 8.58197 0.46835 7.75068
H 13.34787 -1.01534 3.79796
H 14.45697 -1.37309 2.51690
H 13.37823 -2.67627 3.06898
H 14.84181 -2.79638 0.84312
H 14.52169 -4.08722 -0.45224
H 13.54922 -3.96140 1.04171
H 9.98872 -0.41027 -3.99842
H 8.62895 2.23068 -4.18057
H 11.55667 0.83977 -5.05107
H 11.76091 1.74393 -3.50938
H 10.68379 2.32303 -4.67560
H 8.01884 0.22309 -5.58383
H 6.65873 1.06970 -4.87086
H 7.93376 -1.49592 -4.25930
H 6.17082 -1.15048 -4.39174
H 7.11739 -0.72415 -2.89623
H 3.51889 3.83480 -0.84423
H 4.66178 3.65975 -2.21513
H 4.90085 5.01315 -1.03894
H 6.07205 4.28626 3.90115
H 3.40407 3.55223 7.20658
H 4.79747 2.34148 7.54784
H 3.36520 2.00475 6.54413

