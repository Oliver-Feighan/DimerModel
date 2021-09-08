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
Mg 13.13210 1.18897 1.00379
C 15.30118 2.11396 3.70908
C 11.07125 -0.21235 3.33270
C 11.43063 -0.03058 -1.38056
C 15.84610 1.62648 -1.13830
N 13.27577 1.03540 3.20561
C 14.11470 1.47799 4.13009
C 13.61792 1.29618 5.58945
C 12.60262 0.07572 5.34756
C 12.26208 0.32821 3.88026
C 13.26412 -1.28805 5.45653
C 12.87118 2.56622 6.12925
C 12.63153 2.57745 7.64060
H 12.69916 3.85397 8.34659
N 11.27669 0.48514 0.97881
C 10.58582 -0.02500 2.01984
C 9.27669 -0.42139 1.47835
C 9.31087 -0.25543 -0.00016
C 10.71145 0.05993 -0.19771
C 8.12396 -0.82591 2.35071
C 8.25999 -0.37497 -1.06203
O 8.56268 -0.42237 -2.30661
C 6.87350 -0.52886 -0.59787
N 13.62106 0.64126 -0.94127
C 12.68114 0.18659 -1.78847
C 13.20768 0.11912 -3.20218
C 14.62715 0.85499 -3.16753
C 14.71975 1.15510 -1.65533
C 13.37652 -1.28469 -3.78288
C 14.83807 2.02470 -4.16751
C 13.94248 3.18580 -3.79991
N 15.15041 1.93877 1.18357
C 16.12776 2.00962 0.16679
C 17.37781 2.48196 0.79663
C 17.04452 2.63232 2.15118
C 15.68847 2.25937 2.35497
C 18.62415 2.66796 0.13471
C 17.53448 3.01160 3.45686
O 18.58095 3.49697 3.88471
C 16.46327 2.57221 4.54583
C 16.24237 3.73404 5.47527
O 15.72459 4.78585 5.06029
O 16.47629 3.33757 6.73957
C 16.25960 4.45251 7.69967
H 10.47063 -0.87825 3.95294
H 10.85832 -0.34116 -2.25776
H 16.70985 1.88198 -1.75581
H 14.38502 1.00690 6.29772
H 11.71473 0.17185 5.96706
H 14.33588 -1.20868 5.56225
H 13.01172 -1.81122 4.54349
H 12.81046 -1.86325 6.26873
H 11.83264 2.54875 5.77446
H 13.38400 3.46223 5.83803
H 13.23151 1.81948 8.15873
H 11.59820 2.27703 7.83034
H 8.44647 -1.01266 3.38246
H 7.77938 -1.80218 1.99292
H 7.29267 -0.11915 2.30479
H 6.79829 -1.41532 0.02741
H 6.25928 -0.50835 -1.50019
H 6.62032 0.32701 0.03774
H 12.55381 0.64880 -3.88408
H 15.45019 0.16391 -3.40173
H 12.81394 -1.43870 -4.72301
H 13.06981 -2.00295 -3.03417
H 14.40496 -1.58202 -3.98066
H 14.54987 1.71754 -5.17585
H 15.86520 2.37041 -4.21739
H 12.98725 2.91602 -4.26017
H 14.33094 4.11478 -4.21239
H 13.78448 3.23325 -2.71620
H 19.02969 3.60513 0.49213
H 18.50324 2.67280 -0.94655
H 19.16408 1.78357 0.48342
H 16.86377 1.70166 5.06670
H 17.18361 4.56218 8.27575
H 15.34389 4.18909 8.26166
H 16.14500 5.44999 7.29181

