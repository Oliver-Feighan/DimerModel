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
Mg -2.33100 34.27800 27.06700
C  -3.29300 32.73500 29.99200
C  -0.77900 36.82100 29.04000
C  -1.68900 35.99300 24.34500
C  -4.00500 31.86100 25.19100
N  -2.10600 34.64700 29.26400
C  -2.54800 33.87100 30.25800
C  -2.04400 34.45500 31.61800
C  -1.35700 35.84600 31.23800
C  -1.46000 35.80800 29.70100
C  -2.07400 37.11000 31.67200
C  -1.16200 33.50400 32.44100
C  -1.65900 33.38600 33.87200
H  -0.69900 33.03500 35.01800
N  -1.29600 36.10900 26.75300
C  -0.70200 36.96400 27.63600
C  -0.06500 38.02800 26.99800
C  -0.40300 37.90700 25.60500
C  -1.15900 36.61100 25.50300
C  0.74500 39.09500 27.72000
C  -0.12100 38.87200 24.46500
O  -0.65300 38.65500 23.35900
C  0.65800 40.10800 24.66200
N  -2.76200 33.93100 25.11000
C  -2.33900 34.76200 24.09300
C  -2.78500 34.31300 22.65200
C  -3.27700 32.84000 22.98500
C  -3.37800 32.87200 24.52100
C  -3.86900 35.12000 21.83400
C  -2.50400 31.61000 22.32400
C  -1.18700 31.25300 22.96100
N  -3.52800 32.72900 27.43000
C  -4.16000 31.81100 26.61300
C  -4.85300 30.83200 27.34600
C  -4.49800 31.15000 28.66200
C  -3.68000 32.31700 28.68800
C  -5.57400 29.72900 26.69600
C  -4.80700 30.73800 29.97700
O  -5.59800 29.85900 30.32500
C  -4.00900 31.74000 30.96400
C  -4.95800 32.38400 31.85500
O  -5.73200 33.27000 31.57300
O  -4.95600 31.78800 33.04100
C  -5.91200 32.19500 34.02600
H  -0.35300 37.61900 29.65100
H  -1.55300 36.44800 23.36200
H  -4.51800 31.12700 24.56600
H  -2.99900 34.65700 32.10200
H  -0.41300 35.85600 31.78400
H  -2.94800 36.94200 32.30200
H  -2.45800 37.64200 30.80200
H  -1.35600 37.81800 32.08700
H  -0.10200 33.75800 32.45800
H  -1.18900 32.52900 31.95500
H  -2.34600 32.53900 33.87800
H  -2.14000 34.33700 34.10000
H  0.94800 38.81600 28.75400
H  0.13200 39.99100 27.82000
H  1.62000 39.38400 27.13900
H  0.44000 40.83300 23.87700
H  1.70900 39.83900 24.55200
H  0.43200 40.54800 25.63300
H  -1.93500 34.12600 21.99500
H  -4.26000 32.69700 22.53800
H  -4.17100 35.97400 22.43900
H  -4.79700 34.54700 21.82100
H  -3.62200 35.43900 20.82100
H  -2.25200 31.92600 21.31200
H  -3.11500 30.71700 22.19800
H  -1.06000 30.17800 22.83800
H  -0.96900 31.52500 23.99300
H  -0.42200 31.79500 22.40500
H  -5.20900 29.35700 25.73900
H  -6.66000 29.81200 26.67700
H  -5.43800 28.93700 27.43200
H  -3.30800 31.14200 31.54600
H  -5.54400 33.01100 34.64900
H  -6.05600 31.34600 34.69400
H  -6.83100 32.52000 33.53800


