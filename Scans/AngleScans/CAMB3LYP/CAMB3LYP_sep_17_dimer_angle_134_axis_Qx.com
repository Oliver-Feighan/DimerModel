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
Mg 8.92581 0.80689 0.68460
C 10.47462 2.80315 3.23317
C 7.88497 -1.28770 3.16732
C 7.94006 -1.20513 -1.56206
C 11.01052 2.38229 -1.62144
N 9.20578 0.85815 2.87805
C 9.75939 1.69803 3.73972
C 9.46577 1.37712 5.22982
C 9.15860 -0.19215 5.08184
C 8.68535 -0.21272 3.62985
C 10.40022 -1.06142 5.19187
C 8.21848 2.15768 5.77481
C 8.05840 2.13217 7.29648
H 7.52639 3.31905 7.96037
N 7.64201 -0.70309 0.78643
C 7.32206 -1.42682 1.87965
C 6.34876 -2.43225 1.42583
C 6.24469 -2.34883 -0.05637
C 7.31025 -1.40940 -0.34298
C 5.56762 -3.29472 2.37420
C 5.34458 -3.01527 -1.05229
O 5.58699 -2.97696 -2.31031
C 4.22301 -3.79252 -0.50508
N 9.54745 0.45999 -1.26875
C 8.91415 -0.43497 -2.04716
C 9.35592 -0.31545 -3.48624
C 10.24325 1.01330 -3.55431
C 10.23424 1.40049 -2.05919
C 10.16085 -1.49330 -4.03502
C 9.82615 2.08540 -4.59801
C 8.49504 2.68991 -4.21319
N 10.33552 2.44366 0.72658
C 11.11923 2.92214 -0.34624
C 12.00736 3.97028 0.19695
C 11.69243 4.01315 1.56351
C 10.69387 3.04509 1.85544
C 12.98363 4.69772 -0.54050
C 11.98524 4.64988 2.82745
O 12.68164 5.60077 3.18038
C 11.30014 3.80759 3.98830
C 10.57929 4.76711 4.89507
O 9.60301 5.41635 4.48039
O 11.02148 4.59964 6.15474
C 10.32809 5.52177 7.09320
H 7.70392 -2.12657 3.83993
H 7.55761 -1.79870 -2.39566
H 11.62016 2.98880 -2.29462
H 10.30234 1.53059 5.90092
H 8.35801 -0.50248 5.74825
H 11.30322 -0.47066 5.23357
H 10.39903 -1.68897 4.31027
H 10.31087 -1.74004 6.04499
H 7.30557 1.62390 5.48087
H 8.22350 3.17306 5.42896
H 8.96834 1.78504 7.80091
H 7.30658 1.38212 7.55378
H 5.97751 -3.24836 3.39068
H 5.72495 -4.33373 2.06502
H 4.49726 -3.07872 2.35574
H 4.60832 -4.57130 0.14875
H 3.64299 -4.11788 -1.37087
H 3.61108 -3.13181 0.11912
H 8.50320 -0.20268 -4.14445
H 11.28863 0.79224 -3.81538
H 9.70896 -1.94852 -4.93649
H 10.26686 -2.22976 -3.24958
H 11.19694 -1.26907 -4.28286
H 9.68580 1.62460 -5.57888
H 10.55575 2.87951 -4.71669
H 7.77308 1.96984 -4.60990
H 8.37093 3.66788 -4.67386
H 8.37332 2.71274 -3.12401
H 12.89851 5.73182 -0.23407
H 12.83620 4.58657 -1.61273
H 13.89597 4.20211 -0.19804
H 12.09009 3.26611 4.51017
H 11.10431 6.09284 7.61190
H 9.67481 4.88040 7.71417
H 9.73100 6.31786 6.66405

