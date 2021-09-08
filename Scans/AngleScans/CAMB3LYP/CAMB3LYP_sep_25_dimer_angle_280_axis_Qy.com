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
Mg 13.13212 1.18676 1.00550
C 10.76028 -1.48720 1.32700
C 10.52231 3.19511 0.11668
C 15.24358 3.44083 0.27863
C 15.62732 -1.23732 0.79684
N 10.96722 0.78973 0.78680
C 10.17212 -0.24887 0.99528
C 8.65617 0.07593 0.91812
C 8.73726 1.38271 -0.01174
C 10.14973 1.84546 0.33963
C 8.69370 1.06760 -1.49773
C 8.04380 0.40843 2.32399
C 6.51440 0.43202 2.36998
H 5.84705 -0.04129 3.57957
N 12.89462 3.14404 0.77925
C 11.76111 3.79213 0.43806
C 12.11436 5.21917 0.38379
C 13.58678 5.34493 0.55947
C 13.97657 3.94961 0.52436
C 11.08773 6.30747 0.26103
C 14.49533 6.52625 0.71831
O 15.76825 6.41506 0.61759
C 13.84619 7.83249 0.90228
N 15.11509 1.11393 0.38512
C 15.81983 2.24203 0.18850
C 17.29007 1.94366 0.01602
C 17.45992 0.39572 0.38021
C 15.97958 0.02934 0.62462
C 17.86281 2.20429 -1.37697
C 18.50011 0.04101 1.47788
C 18.03637 0.56635 2.81742
N 13.23899 -0.95905 1.23388
C 14.37910 -1.77765 1.07855
C 13.93173 -3.17851 1.21956
C 12.54759 -3.07803 1.42586
C 12.15658 -1.71216 1.39412
C 14.75852 -4.33143 1.10423
C 11.32674 -3.81746 1.65268
O 11.05235 -4.99884 1.85948
C 10.09584 -2.82806 1.47271
C 9.16622 -3.02636 2.63847
O 9.52666 -2.73142 3.79159
O 7.93793 -3.32406 2.17689
C 6.97766 -3.51902 3.29569
H 9.81502 3.85151 -0.39103
H 16.03009 4.18471 0.13190
H 16.35991 -2.04598 0.84211
H 8.05241 -0.68498 0.43826
H 8.00580 2.12800 0.28971
H 8.73454 0.00497 -1.68560
H 9.55491 1.56134 -1.92851
H 7.81762 1.53477 -1.95661
H 8.25530 1.45755 2.56762
H 8.41745 -0.26448 3.07112
H 6.06811 -0.02783 1.47986
H 6.18201 1.47216 2.32953
H 10.10554 5.90460 -0.01537
H 11.37828 6.92605 -0.59501
H 11.03415 6.93829 1.15090
H 13.20060 8.03792 0.05171
H 14.65792 8.53875 1.08713
H 13.19783 7.77987 1.78402
H 17.88726 2.53333 0.70067
H 17.79021 -0.19327 -0.48802
H 18.71580 2.90890 -1.37336
H 17.06683 2.57406 -2.00960
H 18.19173 1.31660 -1.91457
H 19.45482 0.52898 1.26676
H 18.69390 -1.02315 1.56156
H 18.35902 1.61162 2.79769
H 18.51403 0.02108 3.62894
H 16.94228 0.56031 2.88518
H 14.47585 -4.99888 1.90741
H 15.81361 -4.07079 1.15559
H 14.46973 -4.67977 0.10898
H 9.61797 -3.06915 0.52240
H 6.53317 -4.51064 3.16607
H 6.29296 -2.65105 3.25638
H 7.38431 -3.59877 4.29705

