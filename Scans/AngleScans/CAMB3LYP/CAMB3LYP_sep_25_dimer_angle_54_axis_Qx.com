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
Mg 13.13857 1.19155 1.00409
C 14.93544 0.02096 3.88160
C 10.51322 1.82610 3.08722
C 11.37801 1.80444 -1.56344
C 15.41871 -0.56436 -0.96148
N 12.87920 0.91738 3.18349
C 13.65209 0.52205 4.18385
C 13.04398 0.74708 5.59421
C 11.49086 0.76514 5.18634
C 11.62384 1.23221 3.73829
C 10.85247 -0.61392 5.18667
C 13.48799 2.11320 6.22555
C 13.18522 2.26036 7.71834
H 14.13227 2.99629 8.55139
N 11.42675 2.17987 0.82675
C 10.48819 2.35922 1.77975
C 9.41255 3.13575 1.14400
C 9.71763 3.25872 -0.30739
C 10.86710 2.38225 -0.41032
C 8.28037 3.74401 1.91962
C 9.07277 4.02675 -1.42116
O 9.35954 3.79998 -2.64965
C 8.02787 4.98675 -1.03637
N 13.23366 0.52457 -0.96275
C 12.37939 0.98797 -1.89195
C 12.81038 0.58054 -3.28070
C 14.26971 -0.05088 -3.10973
C 14.39873 0.02264 -1.57241
C 11.90221 -0.42258 -3.99160
C 15.40580 0.55758 -3.97676
C 15.69023 1.97481 -3.53410
N 14.97242 0.10186 1.34710
C 15.75298 -0.57810 0.38665
C 16.84010 -1.26346 1.11526
C 16.60129 -0.95105 2.46205
C 15.43519 -0.14567 2.56736
C 17.84178 -2.08889 0.53095
C 17.06441 -1.12823 3.81956
O 18.05397 -1.64196 4.33973
C 15.93423 -0.61154 4.81057
C 16.58705 0.26798 5.84149
O 17.10954 1.34735 5.51196
O 16.29759 -0.20192 7.06859
C 16.91389 0.64370 8.12543
H 9.55680 1.84953 3.61030
H 10.87016 2.07602 -2.49178
H 16.22675 -1.05140 -1.51152
H 13.23305 -0.05027 6.30291
H 10.93696 1.49371 5.77264
H 11.58076 -1.39497 5.34696
H 10.38796 -0.72328 4.21543
H 10.03839 -0.65228 5.91614
H 12.85451 2.91527 5.82540
H 14.53049 2.29367 6.04878
H 12.92917 1.30069 8.18356
H 12.28056 2.86275 7.83120
H 8.23238 3.34556 2.94057
H 7.35140 3.40204 1.45064
H 8.30081 4.83586 1.90841
H 7.23447 4.46418 -0.50731
H 7.75000 5.50136 -1.95824
H 8.45668 5.70602 -0.32970
H 12.87455 1.44201 -3.93404
H 14.28567 -1.11803 -3.37574
H 11.52722 -0.05739 -4.96633
H 11.07818 -0.66330 -3.33315
H 12.34628 -1.39992 -4.17277
H 15.09432 0.61569 -5.02259
H 16.32649 -0.01548 -3.94632
H 14.92717 2.55569 -4.06079
H 16.69208 2.27591 -3.83357
H 15.51209 2.09461 -2.45913
H 18.77940 -1.81880 0.99831
H 17.88316 -1.96008 -0.54862
H 17.46868 -3.07724 0.81242
H 15.46526 -1.48851 5.25841
H 17.52217 -0.01802 8.74977
H 16.07275 1.16754 8.61720
H 17.64925 1.37684 7.81529

