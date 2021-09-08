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
Mg 8.94169 0.80588 0.68604
C 7.65728 -1.18413 3.38219
C 8.74473 3.49503 2.77512
C 9.64086 2.65791 -1.79359
C 7.96968 -1.74756 -1.47750
N 8.27214 1.00315 2.78486
C 7.86349 0.16768 3.72791
C 7.74157 0.80127 5.13974
C 7.54715 2.33342 4.69975
C 8.26132 2.29271 3.35032
C 6.09480 2.71115 4.46031
C 9.04300 0.61000 5.99514
C 8.88684 0.92890 7.48366
H 9.63422 0.12791 8.44922
N 9.62220 2.66918 0.62629
C 9.47159 3.62116 1.57092
C 10.13594 4.82213 1.04131
C 10.54776 4.54980 -0.36260
C 9.92073 3.26084 -0.57593
C 10.39569 6.04296 1.87529
C 11.35949 5.32736 -1.35390
O 11.39225 5.00879 -2.59498
C 12.04013 6.52759 -0.84628
N 8.63400 0.59873 -1.35963
C 9.08241 1.52630 -2.22360
C 8.99206 1.03374 -3.64825
C 8.61349 -0.51653 -3.54363
C 8.45611 -0.63578 -2.01184
C 7.97248 1.74758 -4.53553
C 9.54967 -1.52022 -4.27085
C 10.90294 -1.54458 -3.59747
N 8.15762 -1.19642 0.89601
C 7.79548 -2.08377 -0.14109
C 7.20642 -3.27985 0.49518
C 7.24596 -2.99334 1.86807
C 7.80064 -1.70093 2.07180
C 6.68124 -4.41242 -0.18862
C 6.93481 -3.48595 3.19065
O 6.52887 -4.55419 3.64652
C 7.07081 -2.28449 4.22235
C 7.87224 -2.77083 5.39866
O 9.06701 -3.08771 5.26185
O 7.16344 -2.57593 6.52560
C 7.92118 -3.03259 7.72105
H 8.51011 4.43789 3.26995
H 9.96523 3.21003 -2.67870
H 7.73238 -2.62847 -2.07780
H 6.88595 0.46812 5.71480
H 8.05841 3.00807 5.38153
H 5.44161 1.85271 4.51109
H 6.06447 3.15174 3.47238
H 5.79197 3.50207 5.15246
H 9.77289 1.37990 5.71362
H 9.43413 -0.38183 5.87717
H 7.83377 1.00495 7.78064
H 9.29031 1.92737 7.66882
H 9.83331 6.01410 2.81672
H 9.97407 6.89624 1.33308
H 11.46031 6.22088 2.04141
H 11.30415 7.21110 -0.42963
H 12.64048 6.89758 -1.67969
H 12.69903 6.23369 -0.02161
H 9.94627 1.12927 -4.15170
H 7.62477 -0.72474 -3.97828
H 8.41911 2.18619 -5.44779
H 7.48321 2.51185 -3.94644
H 7.13382 1.13432 -4.86076
H 9.72146 -1.19949 -5.30127
H 9.15405 -2.52967 -4.30949
H 11.41574 -0.68689 -4.04298
H 11.42389 -2.47411 -3.81802
H 10.81106 -1.35167 -2.52230
H 7.03496 -5.28725 0.34059
H 6.99048 -4.42575 -1.23167
H 5.60944 -4.22577 -0.08059
H 6.06258 -1.98399 4.51002
H 7.28619 -3.75310 8.24578
H 8.19778 -2.11224 8.26895
H 8.81514 -3.62167 7.55290

