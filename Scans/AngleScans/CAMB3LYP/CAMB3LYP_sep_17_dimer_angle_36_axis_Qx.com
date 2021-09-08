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
Mg 8.93672 0.81248 0.68498
C 10.18415 -0.91871 3.57045
C 6.56628 2.18451 2.71800
C 7.54080 1.99567 -1.90711
C 10.62845 -1.52098 -1.27429
N 8.53082 0.58552 2.84738
C 9.10894 -0.05120 3.85470
C 8.55223 0.32088 5.25512
C 7.09566 0.82807 4.80787
C 7.41600 1.26188 3.37893
C 6.06269 -0.28496 4.74891
C 9.37482 1.46830 5.93980
C 9.08141 1.67000 7.42808
H 10.18043 2.05820 8.30781
N 7.62130 2.28617 0.49402
C 6.75205 2.72705 1.42747
C 5.99152 3.81211 0.78841
C 6.36919 3.86559 -0.65007
C 7.19431 2.67848 -0.75028
C 5.07689 4.72449 1.55289
C 6.03182 4.81933 -1.75585
O 6.27636 4.54126 -2.98292
C 5.32224 6.04738 -1.36871
N 8.88831 0.19120 -1.29846
C 8.25145 0.91632 -2.23469
C 8.58271 0.42526 -3.62376
C 9.76869 -0.63079 -3.43466
C 9.86138 -0.63385 -1.89297
C 7.43389 -0.23158 -4.38872
C 11.06621 -0.38589 -4.25267
C 11.75939 0.86361 -3.75907
N 10.33112 -0.79892 1.04149
C 10.89572 -1.66653 0.08099
C 11.69217 -2.67049 0.81595
C 11.51563 -2.32840 2.16513
C 10.65264 -1.20378 2.26508
C 12.40909 -3.75296 0.23251
C 11.85453 -2.66944 3.52805
O 12.61844 -3.47556 4.05753
C 10.90612 -1.84934 4.50496
C 11.76310 -1.23772 5.57922
O 12.60457 -0.36659 5.29696
O 11.30073 -1.62099 6.78325
C 12.11163 -1.03084 7.88142
H 5.64661 2.49193 3.21657
H 7.17386 2.43108 -2.83944
H 11.26474 -2.22253 -1.81802
H 8.46121 -0.51081 5.94344
H 6.77431 1.67967 5.40191
H 6.50798 -1.25649 4.90374
H 5.62067 -0.22414 3.76303
H 5.25212 -0.08480 5.45548
H 9.03435 2.43558 5.54840
H 10.42747 1.32060 5.79580
H 8.52553 0.82716 7.85653
H 8.40391 2.52050 7.53606
H 4.87316 4.33873 2.55942
H 4.10435 4.69735 1.04961
H 5.43415 5.75624 1.57607
H 4.38855 5.78517 -0.87685
H 5.24877 6.64241 -2.28108
H 5.92789 6.58304 -0.62921
H 8.93233 1.23825 -4.24820
H 9.46315 -1.64440 -3.73311
H 7.22374 0.25266 -5.36105
H 6.55376 -0.21916 -3.75964
H 7.56008 -1.29429 -4.58851
H 10.82394 -0.21171 -5.30392
H 11.76310 -1.21655 -4.21614
H 11.23173 1.66351 -4.28711
H 12.81489 0.84579 -4.02298
H 11.59029 1.00964 -2.68588
H 13.36770 -3.79679 0.73217
H 12.52521 -3.62019 -0.84111
H 11.73938 -4.58284 0.47344
H 10.17403 -2.54722 4.91309
H 12.46386 -1.86171 8.50051
H 11.45725 -0.28281 8.36721
H 13.04775 -0.55519 7.61332

