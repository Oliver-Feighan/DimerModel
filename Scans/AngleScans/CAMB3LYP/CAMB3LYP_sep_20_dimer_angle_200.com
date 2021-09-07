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
Mg 10.51878 0.94510 0.79151
C 10.18594 3.95782 -1.12990
C 10.01116 -0.78163 -2.10596
C 10.35772 -1.71128 2.51925
C 9.86214 2.85064 3.63465
N 10.15739 1.61390 -1.28560
C 10.11331 2.78286 -1.90683
C 10.06484 2.68685 -3.45538
C 9.44249 1.21219 -3.58551
C 9.92602 0.62507 -2.26106
C 7.92290 1.19832 -3.58768
C 11.49121 2.76765 -4.10420
C 11.49735 2.99086 -5.61810
H 12.55623 3.81382 -6.19602
N 10.71471 -0.95375 0.24879
C 10.48950 -1.47788 -0.97426
C 10.75545 -2.91936 -0.84942
C 10.99631 -3.22932 0.58607
C 10.67130 -1.95252 1.18958
C 10.84171 -3.83652 -2.03466
C 11.43298 -4.47316 1.29907
O 11.34712 -4.58433 2.57308
C 11.87656 -5.59164 0.45425
N 9.95444 0.57407 2.75777
C 10.02687 -0.66462 3.27579
C 9.83613 -0.64150 4.77361
C 9.86146 0.90341 5.18614
C 9.97826 1.53813 3.78297
C 8.54955 -1.28495 5.29039
C 10.88573 1.32430 6.27536
C 12.29334 1.18188 5.74272
N 10.28039 3.04623 1.23477
C 9.98437 3.62532 2.48821
C 9.80602 5.07502 2.26682
C 9.98942 5.23669 0.88525
C 10.24279 3.97388 0.28476
C 9.46770 6.03332 3.26346
C 10.01577 6.18703 -0.20326
O 9.95408 7.41291 -0.28785
C 10.01894 5.38420 -1.57506
C 11.08704 5.97304 -2.45538
O 12.28651 5.88484 -2.13900
O 10.53696 6.34865 -3.62444
C 11.55907 6.92876 -4.53581
H 9.63746 -1.41584 -2.91033
H 10.39803 -2.58144 3.17855
H 9.75234 3.52679 4.48521
H 9.41459 3.40851 -3.93512
H 9.87697 0.67259 -4.42302
H 7.51079 2.17501 -3.38183
H 7.63404 0.49569 -2.81702
H 7.55106 0.78312 -4.52878
H 11.96204 1.77722 -4.05839
H 12.08542 3.52149 -3.62538
H 10.51640 3.30859 -5.99189
H 11.67342 2.03097 -6.10967
H 10.45665 -3.35712 -2.94312
H 10.14849 -4.66421 -1.84962
H 11.84509 -4.24197 -2.18135
H 11.07244 -5.87562 -0.22045
H 12.23305 -6.35709 1.14626
H 12.70707 -5.24853 -0.17276
H 10.64899 -1.15245 5.27490
H 8.89618 1.23060 5.59968
H 8.72999 -2.09471 6.02237
H 7.98936 -1.65654 4.44265
H 7.84226 -0.59757 5.75131
H 10.81265 0.66076 7.14054
H 10.74268 2.33813 6.63408
H 12.50952 0.12055 5.89739
H 12.98100 1.81072 6.30461
H 12.32681 1.37124 4.66352
H 10.09898 6.89536 3.09354
H 9.60446 5.63064 4.26491
H 8.41279 6.19638 3.02712
H 9.02568 5.48288 -2.01469
H 11.20527 7.92460 -4.82008
H 11.68890 6.18938 -5.34835
H 12.53068 7.16139 -4.11594

