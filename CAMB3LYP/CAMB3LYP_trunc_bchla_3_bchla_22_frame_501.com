%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 1.29600 7.71300 26.88100
C  1.46900 9.82700 29.67000
C  2.42600 5.06900 28.94700
C  1.47300 5.72500 24.13400
C  0.92500 10.50000 24.87700
N  1.69900 7.49100 29.04900
C  1.51500 8.49600 30.00300
C  1.43800 7.84300 31.43000
C  1.83300 6.38800 31.15200
C  2.05200 6.30500 29.60100
C  3.12000 6.01700 31.93700
C  0.04000 8.06800 31.99300
C  -0.12000 8.13500 33.50300
H  1.05800 7.66300 34.38500
N  1.84800 5.68100 26.53500
C  2.31100 4.77600 27.50900
C  2.63500 3.49600 26.90400
C  2.31300 3.59100 25.49300
C  1.87500 5.04400 25.33500
C  3.03900 2.24400 27.69300
C  2.37500 2.55400 24.38800
O  2.03400 2.84400 23.24900
C  2.64700 1.14500 24.75200
N  1.13700 8.09300 24.78900
C  1.15400 7.08500 23.87300
C  0.77800 7.64200 22.47500
C  0.53200 9.18300 22.70900
C  0.76500 9.29400 24.23800
C  1.89300 7.33600 21.45900
C  -0.89600 9.67800 22.35700
C  -0.86300 10.38200 21.00400
N  1.30900 9.79300 27.10700
C  1.22300 10.82400 26.23200
C  1.50800 12.06900 26.98600
C  1.50300 11.75700 28.31400
C  1.42300 10.34100 28.32600
C  1.55300 13.42200 26.29600
C  1.55800 12.26300 29.69200
O  1.68500 13.41500 30.08900
C  1.44800 11.03500 30.57500
C  2.56000 11.13900 31.53300
O  3.76600 10.89700 31.38900
O  2.04800 11.65700 32.70500
C  2.99200 11.90900 33.77200
H  2.69600 4.26200 29.63100
H  1.45200 4.99400 23.32300
H  0.89300 11.33000 24.16800
H  2.21800 8.40600 31.94400
H  1.06000 5.77700 31.61800
H  4.00400 5.85500 31.32000
H  2.97100 4.99800 32.29300
H  3.38200 6.77100 32.67900
H  -0.56100 7.28300 31.53300
H  -0.37400 8.97400 31.55200
H  -1.02500 7.59600 33.78000
H  -0.31800 9.19800 33.64500
H  3.22100 2.42800 28.75200
H  3.95300 1.82100 27.27400
H  2.28900 1.47200 27.52500
H  1.81900 0.77800 25.35900
H  3.66100 1.02700 25.13400
H  2.61600 0.48900 23.88200
H  -0.10600 7.10300 22.13500
H  1.34500 9.68200 22.18000
H  2.33200 8.19000 20.94400
H  1.31000 6.77100 20.73200
H  2.77000 6.84600 21.88100
H  -1.17200 10.51700 22.99600
H  -1.64500 8.88600 22.34300
H  -1.71200 9.95300 20.47200
H  -0.02300 10.09700 20.37100
H  -0.94100 11.46900 21.03900
H  1.59500 13.38900 25.20700
H  2.38500 14.03800 26.63500
H  0.58900 13.91000 26.44400
H  0.47700 11.06900 31.06800
H  3.81000 12.49700 33.35500
H  3.40600 11.00100 34.21100
H  2.51900 12.68100 34.37800
Mg 9.20000 47.36600 24.88000
C  7.02600 47.76900 27.74700
C  11.57000 48.44700 26.79100
C  10.94600 47.19300 22.09000
C  6.29100 46.96700 22.73800
N  9.27600 47.92100 27.10800
C  8.32000 48.03400 28.09600
C  8.93500 48.52700 29.35600
C  10.25500 49.18700 28.88100
C  10.38000 48.52600 27.48800
C  10.19000 50.72500 28.75600
C  9.27100 47.32000 30.28400
C  8.91000 47.48600 31.81100
H  8.76300 46.21000 32.73200
N  11.12300 47.75000 24.52700
C  11.99300 48.08000 25.44400
C  13.30700 48.11300 24.83000
C  13.09300 47.62600 23.51900
C  11.67000 47.53200 23.33600
C  14.56100 48.44000 25.66100
C  14.15200 47.33300 22.49300
O  13.84400 46.95300 21.40300
C  15.59100 47.42400 22.88400
N  8.59400 47.38800 22.62900
C  9.59400 47.15300 21.78200
C  9.10700 46.72400 20.44900
C  7.56400 46.85800 20.64800
C  7.41500 47.09900 22.12700
C  9.70100 47.44000 19.18300
C  6.80800 45.57800 20.11700
C  6.00900 45.65200 18.78600
N  7.17400 47.35700 25.08600
C  6.10400 47.11400 24.19600
C  4.83100 47.11900 24.94600
C  5.17200 47.34400 26.27300
C  6.54300 47.52700 26.33600
C  3.42100 47.05000 24.31100
C  4.55900 47.48400 27.63300
O  3.43200 47.38700 28.05800
C  5.73500 47.94200 28.63300
C  5.68600 47.23700 29.86100
O  6.32700 46.19300 30.04600
O  4.84600 47.74900 30.84500
C  4.90400 47.16700 32.22100
H  12.46500 48.75400 27.33500
H  11.53100 47.07200 21.17600
H  5.40800 46.83900 22.11000
H  8.32300 49.32300 29.78100
H  11.13800 49.00400 29.49300
H  10.92100 51.19300 29.41600
H  9.18200 51.02300 29.04300
H  10.33700 51.10600 27.74500
H  10.35100 47.19900 30.19800
H  8.78700 46.45500 29.83100
H  8.00100 48.08400 31.87000
H  9.74500 48.08400 32.17800
H  15.18900 47.54900 25.70000
H  14.38200 48.83100 26.66200
H  15.20000 49.13500 25.11700
H  15.85600 46.71400 23.66700
H  15.91500 48.35300 23.35300
H  16.20400 47.15000 22.02500
H  9.34100 45.66500 20.34000
H  7.22100 47.78300 20.18400
H  10.49600 48.09000 19.54700
H  9.00000 47.94700 18.51900
H  10.12400 46.67900 18.52700
H  6.07000 45.27200 20.85900
H  7.53000 44.76900 20.22800
H  6.23700 44.78300 18.16900
H  6.25500 46.57800 18.26600
H  4.93600 45.62800 18.98000
H  3.49200 46.43400 23.41500
H  3.14200 48.04200 23.95500
H  2.68000 46.64400 25.00000
H  5.61100 49.00700 28.82700
H  5.88500 47.32300 32.67200
H  4.47300 46.16700 32.25100
H  4.32200 47.83100 32.86000


