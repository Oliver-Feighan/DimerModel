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
Mg 15.75611 1.42204 1.19625
C 14.57974 -1.39853 -0.68513
C 14.85313 3.37482 -1.45061
C 17.18900 3.86067 2.63407
C 17.57209 -0.81733 3.15409
N 14.84329 0.93475 -0.75863
C 14.38324 -0.16385 -1.33811
C 13.60134 0.08333 -2.65610
C 14.25178 1.49481 -3.06048
C 14.64582 1.98959 -1.67030
C 15.51466 1.35856 -3.89469
C 12.05543 0.21007 -2.41918
C 11.20275 0.14314 -3.68812
H 9.89778 -0.50831 -3.61624
N 15.56293 3.37613 0.90855
C 15.16257 3.99667 -0.22102
C 15.20818 5.43726 0.07354
C 15.82539 5.62244 1.41499
C 16.24453 4.26521 1.70194
C 14.62382 6.47383 -0.84160
C 16.01921 6.82467 2.28865
O 16.79321 6.79713 3.30993
C 15.34954 8.06134 1.86009
N 17.33689 1.53431 2.54194
C 17.72733 2.71518 3.05283
C 18.69063 2.52166 4.19957
C 18.67871 0.95300 4.51062
C 17.73700 0.47841 3.38219
C 20.12694 2.97909 3.94624
C 18.35746 0.52722 5.96945
C 16.92064 0.86131 6.30008
N 15.90296 -0.72641 1.37131
C 16.74505 -1.45595 2.23891
C 16.57345 -2.88624 1.91140
C 15.65344 -2.88816 0.85212
C 15.29351 -1.55216 0.52800
C 17.25974 -3.96933 2.52949
C 14.91301 -3.71602 -0.07260
O 14.74925 -4.92708 -0.21492
C 14.28121 -2.78194 -1.19274
C 12.83780 -3.17007 -1.36135
O 12.02441 -2.99626 -0.43680
O 12.61208 -3.47780 -2.65159
C 11.19132 -3.85868 -2.87110
H 14.81719 4.04722 -2.30820
H 17.63143 4.65924 3.23409
H 18.02875 -1.58371 3.78401
H 13.78289 -0.64672 -3.43573
H 13.51362 2.15679 -3.50576
H 15.83306 0.33011 -3.97861
H 16.26742 1.94746 -3.38716
H 15.37314 1.82818 -4.87232
H 11.82592 1.23243 -2.09230
H 11.71606 -0.52595 -1.71653
H 11.77195 -0.23048 -4.54798
H 10.92393 1.16057 -3.97291
H 14.38730 6.05325 -1.82678
H 15.41377 7.20591 -1.04132
H 13.76732 6.98924 -0.40180
H 15.69266 8.33135 0.86419
H 15.53254 8.78419 2.65751
H 14.27406 7.86732 1.78059
H 18.35653 3.05701 5.07983
H 19.65835 0.49179 4.31716
H 20.48384 3.72435 4.68195
H 20.18620 3.37747 2.94209
H 20.86809 2.18201 3.92267
H 18.97659 1.08918 6.67303
H 18.52972 -0.52689 6.15962
H 16.97117 1.91766 6.58057
H 16.56574 0.24984 7.12723
H 16.28446 0.78653 5.41051
H 16.52464 -4.74307 2.70686
H 17.74206 -3.65867 3.45394
H 17.98589 -4.21107 1.74876
H 14.85538 -2.93325 -2.10762
H 11.19423 -4.85066 -3.33333
H 10.74728 -3.03158 -3.45631
H 10.57869 -4.03548 -1.99481

