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
Mg 10.51770 0.94462 0.80637
C 10.62788 -2.64179 0.87378
C 9.40091 0.82195 4.02698
C 9.98069 4.06704 0.63431
C 10.53118 0.74118 -2.67273
N 10.10800 -0.75054 2.16681
C 10.24525 -2.06640 2.10340
C 10.03400 -2.80096 3.45445
C 9.11490 -1.70421 4.18332
C 9.59006 -0.45589 3.44265
C 7.63092 -1.88546 3.91074
C 11.38004 -3.02818 4.22813
C 11.29029 -4.00688 5.40109
H 12.43599 -4.86916 5.67734
N 10.27116 2.28667 2.24745
C 9.83080 2.06089 3.50299
C 9.80796 3.37506 4.16374
C 10.11345 4.41710 3.14613
C 10.11312 3.61501 1.93918
C 9.59462 3.53973 5.64049
C 10.35211 5.89362 3.24173
O 10.36463 6.64242 2.20154
C 10.47774 6.45106 4.59634
N 10.07552 2.21395 -0.77971
C 9.94208 3.53823 -0.58888
C 9.90061 4.27755 -1.90496
C 10.27853 3.19380 -3.01863
C 10.39078 1.93998 -2.12380
C 8.56385 4.92776 -2.26088
C 11.46410 3.53755 -3.96146
C 12.75722 3.55690 -3.17862
N 10.75479 -0.62531 -0.65895
C 10.70119 -0.49577 -2.06410
C 10.80208 -1.85628 -2.63090
C 10.88488 -2.69157 -1.50665
C 10.81851 -1.90873 -0.32248
C 10.76056 -2.18604 -4.01497
C 11.00264 -4.05845 -1.05226
O 11.18438 -5.14402 -1.60202
C 10.71421 -4.09817 0.51030
C 11.79299 -4.91995 1.16102
O 12.97351 -4.52881 1.15887
O 11.22438 -5.91712 1.86294
C 12.25367 -6.75357 2.53580
H 8.83247 0.88908 4.95501
H 9.90420 5.15039 0.51608
H 10.64126 0.60181 -3.75029
H 9.50256 -3.74226 3.38149
H 9.35028 -1.63217 5.24199
H 7.44748 -2.65337 3.17398
H 7.27814 -0.92671 3.55383
H 7.09616 -2.07424 4.84600
H 11.64139 -2.10974 4.76934
H 12.15739 -3.33848 3.55726
H 10.36414 -4.59397 5.37812
H 11.22046 -3.43229 6.32790
H 9.23227 2.61209 6.10048
H 8.76765 4.24680 5.76722
H 10.47619 3.93156 6.15232
H 9.57328 6.23634 5.16052
H 10.73376 7.50314 4.45681
H 11.29855 5.93888 5.11069
H 10.63546 5.07297 -1.92501
H 9.44303 3.01225 -3.71059
H 8.64301 6.01587 -2.44481
H 7.86313 4.72726 -1.46131
H 8.05758 4.49985 -3.12447
H 11.33840 4.53999 -4.37796
H 11.56470 2.85219 -4.79642
H 12.76610 4.55988 -2.74153
H 13.60793 3.40720 -3.84039
H 12.72775 2.83739 -2.35213
H 11.53557 -2.92197 -4.18307
H 10.90522 -1.30581 -4.63791
H 9.74365 -2.58020 -4.09114
H 9.72530 -4.53653 0.64986
H 12.08523 -7.78643 2.21565
H 12.15377 -6.53848 3.61632
H 13.28714 -6.60760 2.24412

