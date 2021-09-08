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
Mg 13.12983 1.18513 0.99250
C 13.02139 -1.59148 -1.27856
C 13.33129 3.21871 -1.73853
C 13.67954 3.61070 2.96274
C 14.06185 -1.06732 3.48319
N 13.16250 0.75130 -1.17614
C 13.05343 -0.33719 -1.92289
C 12.89338 -0.06156 -3.44209
C 13.57534 1.39133 -3.49604
C 13.31347 1.83777 -2.05910
C 15.07839 1.34290 -3.71475
C 11.38893 -0.00672 -3.88438
C 11.16135 -0.04965 -5.39703
H 9.98750 -0.75462 -5.90456
N 12.97012 3.14112 0.69724
C 13.05424 3.79593 -0.47980
C 12.89110 5.22237 -0.15919
C 12.86859 5.37077 1.32142
C 13.19999 4.01952 1.72684
C 12.69497 6.27549 -1.21080
C 12.60641 6.53889 2.22306
O 12.87331 6.49563 3.47595
C 12.11544 7.76635 1.57981
N 13.98019 1.29854 2.88627
C 14.05098 2.46953 3.54332
C 14.44474 2.26237 4.98638
C 14.38802 0.68198 5.22457
C 14.04271 0.22225 3.79111
C 15.82503 2.78707 5.38120
C 13.50128 0.17821 6.39613
C 12.04417 0.44009 6.09013
N 13.30623 -0.96159 1.16146
C 13.73855 -1.69592 2.28748
C 13.80128 -3.11555 1.88357
C 13.42021 -3.10574 0.53332
C 13.15928 -1.77166 0.11903
C 14.21829 -4.19764 2.70898
C 13.18963 -3.91940 -0.63864
O 13.16870 -5.12828 -0.86647
C 13.04339 -2.96117 -1.89832
C 11.83215 -3.39785 -2.67597
O 10.69480 -3.29830 -2.18290
O 12.19337 -3.65547 -3.94622
C 11.02383 -4.08172 -4.76004
H 13.62627 3.92734 -2.51295
H 13.78044 4.39808 3.71329
H 14.24893 -1.84285 4.22899
H 13.42884 -0.74782 -4.08715
H 13.06119 2.04319 -4.19754
H 15.45828 0.33266 -3.67965
H 15.51053 1.93742 -2.92049
H 15.34012 1.85039 -4.64764
H 10.98645 0.98993 -3.66204
H 10.82417 -0.78663 -3.41171
H 12.06157 -0.36121 -5.94062
H 10.97451 0.96785 -5.74886
H 12.92300 5.89099 -2.21258
H 13.45337 7.04644 -1.03664
H 11.70586 6.73629 -1.16614
H 12.83380 8.09451 0.83233
H 11.90222 8.45898 2.39631
H 11.18813 7.53388 1.04439
H 13.73939 2.74371 5.65263
H 15.38074 0.26899 5.45647
H 15.79393 3.51173 6.21665
H 16.28337 3.23256 4.50827
H 16.54854 2.02171 5.65682
H 13.73084 0.73169 7.30998
H 13.63406 -0.87582 6.61607
H 11.91257 1.48380 6.39083
H 11.40564 -0.22175 6.67176
H 11.85145 0.38067 5.01264
H 13.52125 -5.00700 2.53701
H 14.24420 -3.91041 3.75808
H 15.21945 -4.37494 2.30709
H 13.95925 -3.04812 -2.48415
H 11.27740 -5.05085 -5.20069
H 10.82582 -3.24710 -5.45857
H 10.10765 -4.32205 -4.23343

