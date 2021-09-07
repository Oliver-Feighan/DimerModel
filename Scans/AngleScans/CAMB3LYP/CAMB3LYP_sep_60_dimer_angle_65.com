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
Mg 31.53962 2.83539 2.40245
C 31.70918 -0.65140 1.57020
C 30.50536 1.88967 5.51220
C 30.94830 5.88969 3.02612
C 31.46990 3.50572 -1.01683
N 31.19121 0.84616 3.30458
C 31.34801 -0.40907 2.91194
C 31.18227 -1.46158 4.04080
C 30.26403 -0.60136 5.03867
C 30.70047 0.80220 4.62369
C 28.77685 -0.74126 4.75928
C 32.55060 -1.84504 4.70610
C 32.50581 -3.08689 5.59916
H 33.67176 -3.96562 5.62837
N 31.30740 3.77026 4.13758
C 30.90208 3.22917 5.30557
C 30.87451 4.33639 6.27380
C 31.13780 5.60556 5.54263
C 31.12039 5.12975 4.17391
C 30.69532 4.12317 7.74887
C 31.35496 7.01651 5.99899
O 31.32956 8.00099 5.17863
C 31.50523 7.22139 7.44707
N 31.03778 4.44998 1.19267
C 30.88781 5.68169 1.71068
C 30.80175 6.72457 0.62193
C 31.16927 5.96105 -0.73437
C 31.32390 4.52654 -0.18331
C 29.44626 7.41368 0.46648
C 32.32535 6.55467 -1.58511
C 33.63705 6.40650 -0.84832
N 31.76547 1.68571 0.58716
C 31.67488 2.16013 -0.73979
C 31.78354 0.98632 -1.63019
C 31.90770 -0.10077 -0.75195
C 31.85819 0.36062 0.59124
C 31.71295 1.01107 -3.05164
C 32.05870 -1.53486 -0.65564
O 32.24418 -2.44496 -1.46252
C 31.80988 -1.96901 0.85300
C 32.91760 -2.90331 1.25620
O 34.09125 -2.49835 1.32809
O 32.38274 -4.05609 1.69822
C 33.44177 -5.01120 2.12019
H 29.95915 1.71100 6.43881
H 30.85146 6.96640 3.18364
H 31.55539 3.64172 -2.09705
H 30.66441 -2.36634 3.74584
H 30.52446 -0.79031 6.07685
H 28.58753 -1.30515 3.85798
H 28.39991 0.26826 4.66016
H 28.26860 -1.16877 5.62833
H 32.81050 -1.08498 5.45410
H 33.31594 -1.96135 3.96356
H 31.58895 -3.66977 5.44893
H 32.44980 -2.76305 6.64128
H 30.35950 3.10249 7.96989
H 29.86049 4.75815 8.06465
H 31.58292 4.39420 8.32453
H 30.61864 6.85318 7.95780
H 31.74073 8.28036 7.56950
H 32.34671 6.61520 7.80071
H 31.52297 7.51571 0.78636
H 30.31986 5.93951 -1.43284
H 29.50329 8.51476 0.55846
H 28.76895 7.00502 1.20457
H 28.92564 7.20355 -0.46625
H 32.17321 7.62629 -1.73557
H 32.41619 6.10138 -2.56657
H 33.64064 7.26885 -0.17496
H 34.47334 6.44502 -1.54343
H 33.63972 5.50326 -0.22720
H 32.49529 0.35731 -3.41365
H 31.82789 2.02170 -3.43789
H 30.70094 0.62626 -3.20339
H 30.83193 -2.44976 0.89851
H 33.28208 -5.93511 1.55577
H 33.36531 -5.07440 3.22196
H 34.46519 -4.77466 1.85351

