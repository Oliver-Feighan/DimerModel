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
Mg 13.14637 1.17584 0.99962
C 13.32406 -0.37812 -2.23034
C 12.44357 -1.90706 2.27884
C 12.53994 2.57756 3.78078
C 12.72782 4.22125 -0.64232
N 12.93852 -0.80930 0.04660
C 13.08874 -1.28850 -1.17902
C 13.06051 -2.83774 -1.27184
C 12.22163 -3.14337 0.06298
C 12.57944 -1.89550 0.86763
C 10.71845 -3.15320 -0.15993
C 14.49639 -3.46682 -1.20455
C 14.57049 -4.94379 -1.59801
H 15.75600 -5.41640 -2.30770
N 13.05969 0.41373 2.82999
C 12.78374 -0.86348 3.16741
C 12.82179 -0.90848 4.63734
C 12.98009 0.48138 5.14524
C 12.84215 1.22997 3.91209
C 12.79219 -2.19327 5.41314
C 13.20412 1.04443 6.51600
O 13.07386 2.29608 6.75936
C 13.48892 0.07535 7.58430
N 12.48806 3.08236 1.50376
C 12.35742 3.45095 2.79023
C 12.13925 4.94027 2.91149
C 12.39290 5.53223 1.44765
C 12.63719 4.21420 0.68052
C 10.75778 5.37243 3.40258
C 13.44527 6.66756 1.32025
C 14.82568 6.12744 1.61709
N 13.22446 1.86288 -1.04731
C 12.99330 3.17468 -1.51616
C 13.04520 3.12511 -2.99172
C 13.28197 1.77300 -3.28194
C 13.35123 1.02877 -2.07330
C 12.83634 4.22648 -3.86893
C 13.47850 0.80429 -4.33630
O 13.60814 0.86464 -5.55830
C 13.38856 -0.65260 -3.70720
C 14.55388 -1.45307 -4.22083
O 15.71814 -1.14454 -3.91114
O 14.09431 -2.57185 -4.81032
C 15.21354 -3.40436 -5.32613
H 11.99468 -2.78209 2.74968
H 12.43116 3.12978 4.71710
H 12.70451 5.14518 -1.22412
H 12.53993 -3.23160 -2.13651
H 12.58639 -4.03694 0.56281
H 10.45712 -2.83523 -1.15830
H 10.30780 -2.47194 0.57377
H 10.30815 -4.13818 0.08050
H 14.80801 -3.53628 -0.15440
H 15.18876 -2.90008 -1.79619
H 13.65869 -5.28140 -2.10556
H 14.60758 -5.54445 -0.68595
H 12.50581 -3.04037 4.77768
H 11.97588 -2.10989 6.13878
H 13.72420 -2.38277 5.94989
H 12.66588 -0.63132 7.65967
H 13.70786 0.67367 8.47084
H 14.37609 -0.50350 7.30416
H 12.85241 5.38102 3.59723
H 11.48053 5.97764 1.02453
H 10.79510 6.00450 4.30987
H 10.16607 4.48456 3.58177
H 10.15465 5.90839 2.67170
H 13.25170 7.44814 2.06011
H 13.45204 7.14278 0.34500
H 14.87241 6.15670 2.70981
H 15.58927 6.76231 1.17222
H 14.91188 5.07824 1.31145
H 13.59642 4.16193 -4.63614
H 12.88718 5.17398 -3.33656
H 11.82462 4.00683 -4.22029
H 12.43235 -1.08147 -4.00939
H 15.02366 -3.57060 -6.39103
H 15.24605 -4.29791 -4.67478
H 16.19981 -2.95613 -5.35462

