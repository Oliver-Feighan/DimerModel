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
Mg 13.13072 1.17991 1.00455
C 13.15161 3.89258 3.35403
C 13.77645 -0.90113 3.62884
C 13.56500 -1.16393 -1.08954
C 13.62584 3.53909 -1.51236
N 13.40607 1.56479 3.16518
C 13.33067 2.62769 3.95188
C 13.37376 2.30836 5.47034
C 14.13506 0.89627 5.39760
C 13.72320 0.46875 3.99038
C 15.64873 1.02607 5.43354
C 11.94029 2.15684 6.09019
C 11.89836 2.15171 7.61989
H 10.75926 2.77606 8.28705
N 13.11496 -0.78831 1.25837
C 13.37855 -1.46416 2.39636
C 13.25494 -2.88986 2.05540
C 13.05872 -3.00513 0.58471
C 13.26390 -1.62805 0.18266
C 13.24703 -3.97648 3.09106
C 12.75153 -4.16530 -0.31298
O 12.85967 -4.07788 -1.58714
C 12.41073 -5.43349 0.34812
N 13.74676 1.15932 -0.98100
C 13.79974 0.00981 -1.67651
C 14.00132 0.27250 -3.14986
C 13.83002 1.85208 -3.33201
C 13.63906 2.25800 -1.85412
C 15.34906 -0.16264 -3.72487
C 12.77996 2.33085 -4.37164
C 11.38790 1.97856 -3.89895
N 13.16818 3.33655 0.87974
C 13.41830 4.12031 -0.26797
C 13.45301 5.53151 0.16762
C 13.24189 5.46866 1.55325
C 13.10669 4.11257 1.95606
C 13.70603 6.65455 -0.66957
C 13.11326 6.24046 2.76835
O 13.05487 7.44055 3.03311
C 13.17520 5.24661 4.00723
C 12.04682 5.59481 4.93891
O 10.86458 5.44136 4.58483
O 12.54689 5.84344 6.16304
C 11.46484 6.18280 7.12523
H 14.20242 -1.60927 4.34017
H 13.61550 -1.92669 -1.87001
H 13.67744 4.34109 -2.25182
H 13.94638 3.00925 6.06586
H 13.74720 0.19990 6.13651
H 15.96601 2.05705 5.38282
H 16.01162 0.47588 4.57526
H 16.05031 0.51311 6.31214
H 11.56826 1.14406 5.88865
H 11.28021 2.91364 5.71314
H 12.84031 2.50195 8.05923
H 11.81171 1.11730 7.96120
H 13.57524 -3.60269 4.06879
H 14.01902 -4.69828 2.80312
H 12.28643 -4.49229 3.15306
H 13.23227 -5.73691 0.99268
H 12.13670 -6.11823 -0.45681
H 11.54501 -5.26718 0.99888
H 13.24662 -0.23315 -3.73964
H 14.76281 2.32676 -3.67027
H 15.25491 -0.86843 -4.57167
H 15.93472 -0.60102 -2.92778
H 15.99057 0.64924 -4.06328
H 12.92518 1.81272 -5.32274
H 12.82727 3.39557 -4.57443
H 11.27721 0.93621 -4.21254
H 10.64773 2.61578 -4.37868
H 11.32613 2.00193 -2.80474
H 12.99249 7.41826 -0.39009
H 13.61836 6.39354 -1.72216
H 14.73787 6.87987 -0.38708
H 14.15004 5.37260 4.47987
H 11.71764 7.15450 7.56063
H 11.39980 5.32228 7.81735
H 10.47919 6.38207 6.72119

