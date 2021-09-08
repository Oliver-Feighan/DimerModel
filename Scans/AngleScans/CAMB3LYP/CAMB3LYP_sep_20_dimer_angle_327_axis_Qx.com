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
Mg 10.51836 0.94679 0.80591
C 8.95976 -0.89261 3.46419
C 10.63353 3.61864 2.92312
C 11.47628 2.72071 -1.64420
C 9.26155 -1.44284 -1.39768
N 9.85188 1.20589 2.89898
C 9.32951 0.41911 3.82766
C 9.26989 1.04903 5.24512
C 9.27496 2.59773 4.82043
C 9.99569 2.48095 3.47894
C 7.88469 3.15726 4.56857
C 10.52590 0.68737 6.11325
C 10.39184 1.00865 7.60340
H 11.02030 0.11064 8.56833
N 11.42803 2.71035 0.77527
C 11.38591 3.66434 1.72890
C 12.20248 4.77751 1.22072
C 12.59492 4.46947 -0.18136
C 11.81390 3.27166 -0.41659
C 12.60262 5.94775 2.07151
C 13.51051 5.14865 -1.15430
O 13.51901 4.84075 -2.39846
C 14.32979 6.24880 -0.62520
N 10.21347 0.80013 -1.24540
C 10.78583 1.67251 -2.09351
C 10.65274 1.20928 -3.52463
C 10.08131 -0.28211 -3.44208
C 9.89049 -0.39576 -1.91366
C 9.74231 2.05426 -4.41541
C 10.89342 -1.38823 -4.16985
C 12.22415 -1.58904 -3.48122
N 9.48658 -0.94312 0.98396
C 9.02934 -1.76767 -0.06729
C 8.28670 -2.88650 0.54845
C 8.34420 -2.62079 1.92489
C 9.05401 -1.41037 2.14973
C 7.63240 -3.93733 -0.15423
C 7.95669 -3.08344 3.23809
O 7.41406 -4.09666 3.67703
C 8.22908 -1.91882 4.28492
C 8.94793 -2.51356 5.46478
O 10.09515 -2.97665 5.33815
O 8.25473 -2.24230 6.58561
C 8.93372 -2.80230 7.78444
H 10.51272 4.57857 3.42593
H 11.87875 3.23641 -2.51916
H 8.92334 -2.28099 -2.01069
H 8.37189 0.82035 5.80644
H 9.85801 3.19606 5.51571
H 7.12835 2.38722 4.60201
H 7.92261 3.60789 3.58545
H 7.67462 3.97309 5.26615
H 11.35020 1.36222 5.84896
H 10.79097 -0.34453 5.98849
H 9.35290 1.21345 7.88907
H 10.91500 1.94665 7.80459
H 12.02898 5.98048 3.00601
H 12.29842 6.85255 1.53424
H 13.67892 5.98888 2.25192
H 13.68009 7.01521 -0.20930
H 14.98249 6.54865 -1.44735
H 14.93592 5.86634 0.20362
H 11.61779 1.18915 -4.01592
H 9.07995 -0.36017 -3.89044
H 10.25216 2.44227 -5.31740
H 9.34526 2.86809 -3.82333
H 8.83758 1.55445 -4.75727
H 11.11736 -1.08146 -5.19449
H 10.37478 -2.33955 -4.22456
H 12.84622 -0.79820 -3.91097
H 12.62713 -2.57441 -3.70633
H 12.14337 -1.39673 -2.40505
H 7.86670 -4.85484 0.36901
H 7.95093 -3.97911 -1.19374
H 6.59119 -3.61859 -0.05646
H 7.26291 -1.49690 4.56434
H 8.20663 -3.44248 8.29356
H 9.31655 -1.92944 8.34596
H 9.74878 -3.49732 7.61991

