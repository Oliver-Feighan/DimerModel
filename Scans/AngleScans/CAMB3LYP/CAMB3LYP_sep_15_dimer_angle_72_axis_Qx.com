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
Mg 7.88060 0.71817 0.60467
C 10.02573 0.25002 3.44339
C 5.24279 0.57449 2.76242
C 5.95082 0.67705 -1.91359
C 10.54085 -0.30708 -1.39974
N 7.77555 0.44495 2.79706
C 8.65857 0.33905 3.77860
C 8.04745 0.40872 5.20382
C 6.55516 -0.06677 4.84909
C 6.49927 0.37341 3.38766
C 6.37576 -1.57491 4.90086
C 8.06249 1.86406 5.79014
C 7.76784 1.95662 7.28893
H 8.46149 2.97468 8.07308
N 5.94279 1.12305 0.46492
C 5.01990 1.03299 1.44544
C 3.74047 1.41895 0.83040
C 3.95469 1.58515 -0.63289
C 5.31610 1.10410 -0.75631
C 2.49604 1.67125 1.63113
C 3.07492 2.08136 -1.74013
O 3.38589 1.91637 -2.97257
C 1.79432 2.68286 -1.34044
N 8.12665 0.05262 -1.34951
C 7.14705 0.20027 -2.25863
C 7.64691 -0.09685 -3.65227
C 9.23399 -0.24070 -3.51789
C 9.37370 -0.08329 -1.98782
C 7.07618 -1.35306 -4.30998
C 10.10273 0.66192 -4.43617
C 9.94537 2.11067 -4.03400
N 9.97009 0.25980 0.90796
C 10.89775 -0.17505 -0.06386
C 12.16225 -0.46793 0.64161
C 11.87339 -0.20304 1.98879
C 10.51827 0.20535 2.11668
C 13.35490 -0.96106 0.04117
C 12.40366 -0.18625 3.33312
O 13.51684 -0.35261 3.82999
C 11.19509 -0.01371 4.35096
C 11.56973 1.05602 5.33981
O 11.72331 2.23327 4.96951
O 11.47206 0.55794 6.58590
C 11.82317 1.58496 7.60264
H 4.34009 0.31738 3.31733
H 5.35993 0.74945 -2.82961
H 11.44548 -0.53734 -1.96656
H 8.49260 -0.26880 5.92256
H 5.81815 0.47273 5.43826
H 7.31413 -2.08728 5.05282
H 5.94298 -1.85253 3.94867
H 5.63284 -1.84037 5.65844
H 7.20150 2.41827 5.39481
H 8.99280 2.35233 5.57394
H 7.83393 0.97971 7.78311
H 6.72429 2.25311 7.41977
H 2.60038 1.30932 2.66143
H 1.70694 1.04454 1.20168
H 2.17684 2.71506 1.59582
H 1.21591 1.95732 -0.77344
H 1.34686 3.05757 -2.26308
H 1.99732 3.52099 -0.66446
H 7.42401 0.72160 -4.32575
H 9.57294 -1.25838 -3.76135
H 6.58134 -1.15204 -5.27889
H 6.38466 -1.81615 -3.61868
H 7.79639 -2.15043 -4.48525
H 9.76156 0.58847 -5.47172
H 11.15614 0.40271 -4.42500
H 9.02648 2.41066 -4.54658
H 10.79648 2.69723 -4.37399
H 9.76681 2.20281 -2.95637
H 14.17446 -0.40006 0.47024
H 13.32632 -0.85929 -1.04169
H 13.31386 -2.00711 0.35638
H 11.03274 -0.97837 4.83330
H 12.62254 1.16332 8.21983
H 10.87416 1.83816 8.11158
H 12.28686 2.49950 7.25184

