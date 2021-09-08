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
Mg 8.92718 0.80010 0.68537
C 8.47747 3.50954 2.99530
C 9.99617 -1.07631 3.32562
C 9.75887 -1.45725 -1.38346
C 8.92812 3.16553 -1.87401
N 9.16147 1.26782 2.83447
C 8.90094 2.31142 3.60736
C 9.02866 2.03266 5.12905
C 10.04034 0.78735 5.06130
C 9.69261 0.26543 3.66869
C 11.50307 1.19935 5.06548
C 7.65976 1.62593 5.77918
C 7.64515 1.64000 7.30932
H 6.42038 2.05129 7.98994
N 9.28567 -1.13137 0.96732
C 9.69053 -1.72561 2.10935
C 9.83125 -3.15498 1.79113
C 9.63557 -3.33092 0.32638
C 9.57166 -1.94710 -0.09912
C 10.04494 -4.20538 2.84205
C 9.53679 -4.54376 -0.54852
O 9.60524 -4.46009 -1.82569
C 9.45140 -5.84150 0.13708
N 9.50278 0.86042 -1.31147
C 9.75910 -0.27084 -1.99154
C 9.88305 -0.00102 -3.47215
C 9.41505 1.51493 -3.67325
C 9.17602 1.90383 -2.19778
C 11.27871 -0.18574 -4.06745
C 8.27650 1.76987 -4.69861
C 6.98350 1.17114 -4.19357
N 8.55676 2.92295 0.52928
C 8.63595 3.71936 -0.63413
C 8.41225 5.11948 -0.21935
C 8.23993 5.04259 1.17093
C 8.36863 3.69252 1.59551
C 8.43575 6.25511 -1.07721
C 7.98897 5.79787 2.37727
O 7.71064 6.97019 2.62610
C 8.25724 4.85527 3.62865
C 7.09929 5.00203 4.57741
O 5.96116 4.62336 4.24896
O 7.56419 5.36158 5.78777
C 6.45391 5.50888 6.76624
H 10.55942 -1.67934 4.03836
H 9.93868 -2.21064 -2.15389
H 8.81576 3.94982 -2.62564
H 9.46932 2.83889 5.70314
H 9.80261 0.04373 5.81752
H 11.62014 2.27048 4.99390
H 11.94842 0.71194 4.20808
H 12.00851 0.78638 5.94316
H 7.48127 0.55795 5.59940
H 6.86304 2.23874 5.40455
H 8.51177 2.16837 7.72496
H 7.76008 0.61387 7.66690
H 10.31341 -3.75949 3.80771
H 10.93383 -4.77455 2.54916
H 9.19949 -4.89102 2.93031
H 10.32599 -5.97405 0.76953
H 9.29743 -6.57957 -0.65249
H 8.58087 -5.82907 0.80240
H 9.22699 -0.64956 -4.03968
H 10.23627 2.15010 -4.03656
H 11.30465 -0.91149 -4.90214
H 11.94956 -0.49238 -3.27595
H 11.75054 0.72595 -4.42994
H 8.50051 1.27148 -5.64496
H 8.11956 2.82083 -4.91737
H 7.06535 0.12116 -4.49011
H 6.12888 1.64968 -4.66756
H 6.93677 1.20180 -3.09879
H 7.59624 6.87625 -0.79454
H 8.38106 5.96378 -2.12407
H 9.41149 6.67489 -0.81839
H 9.19883 5.17017 4.08010
H 6.52694 6.51828 7.18275
H 6.56328 4.66378 7.47164
H 5.44176 5.51261 6.37897

