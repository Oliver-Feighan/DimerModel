%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 1.39700 7.67200 26.33400
C  1.72400 9.74100 29.05100
C  2.01300 4.99100 28.38100
C  1.18200 5.57500 23.58000
C  1.22200 10.43500 24.24700
N  1.88100 7.36900 28.50400
C  1.86000 8.39700 29.40300
C  2.08800 7.91400 30.75300
C  2.40000 6.43300 30.60700
C  2.12200 6.20500 29.05700
C  3.81300 5.91300 31.04400
C  0.86900 8.21500 31.69000
C  1.18400 8.50600 33.19800
H  2.44800 7.85400 33.79200
N  1.46600 5.54900 26.07500
C  1.71900 4.62900 27.07800
C  1.70700 3.27200 26.49500
C  1.45400 3.42700 25.06400
C  1.33100 4.87900 24.84700
C  1.96100 1.99300 27.27600
C  1.22600 2.33500 24.04500
O  0.94200 2.65000 22.87600
C  1.34900 0.88300 24.32000
N  1.30400 7.94800 24.20100
C  1.17000 6.94200 23.33200
C  0.89900 7.47700 21.91600
C  0.64800 9.06900 22.12400
C  1.00600 9.14600 23.65000
C  2.06300 7.12000 20.94100
C  -0.80600 9.55500 21.87400
C  -1.00500 10.45600 20.64600
N  1.52400 9.70800 26.54100
C  1.45100 10.73400 25.59600
C  1.56500 12.05000 26.26300
C  1.69000 11.65900 27.61700
C  1.66900 10.24300 27.72900
C  1.49300 13.37100 25.67200
C  1.73800 12.19000 28.99500
O  1.68200 13.31300 29.42100
C  1.80500 10.95700 29.96500
C  2.97400 11.12000 30.84400
O  4.12200 10.98100 30.39200
O  2.68100 11.28600 32.16300
C  3.77100 11.25200 33.13300
H  2.21100 4.15900 29.06000
H  1.15400 4.90100 22.72100
H  1.03400 11.27000 23.56900
H  2.96600 8.42700 31.14700
H  1.62000 5.81700 31.05600
H  4.56300 6.57700 31.47500
H  4.27400 5.43500 30.18000
H  3.77100 5.07400 31.73900
H  0.12500 7.42500 31.58700
H  0.27300 9.01100 31.24400
H  0.33700 8.21300 33.81900
H  1.41300 9.55400 33.38900
H  2.02800 2.15000 28.35200
H  2.91300 1.54300 26.99300
H  1.08700 1.34400 27.22400
H  1.03600 0.37200 23.41000
H  0.74800 0.54700 25.16600
H  2.36800 0.65900 24.63600
H  -0.06900 7.08000 21.61000
H  1.35900 9.68000 21.56800
H  1.72500 6.67300 20.00600
H  2.74700 6.40500 21.39900
H  2.59200 8.06400 20.80800
H  -1.22600 10.10700 22.71400
H  -1.49800 8.72100 21.75500
H  -0.24900 10.22000 19.89700
H  -0.86300 11.48900 20.96200
H  -2.04500 10.34200 20.34300
H  1.26400 14.14800 26.40100
H  0.76200 13.44400 24.86700
H  2.50500 13.65900 25.38700
H  0.85200 11.05500 30.48600
H  4.47400 12.08100 33.04700
H  4.27000 10.28300 33.12600
H  3.36500 11.36600 34.13800
Mg 16.41100 52.36700 25.56600
C  18.01800 50.88800 28.30000
C  13.92400 53.47100 27.63400
C  14.78500 53.43800 22.89200
C  18.97700 50.98900 23.49400
N  16.01700 52.13900 27.73500
C  16.81400 51.56600 28.70100
C  16.31600 51.93200 30.08700
C  14.87100 52.50500 29.80500
C  15.00600 52.79900 28.27200
C  13.75200 51.51700 30.10900
C  17.38000 52.82700 30.82100
C  16.86600 53.42500 32.16000
H  15.83800 52.61100 32.91800
N  14.67300 53.30500 25.31300
C  13.70900 53.59600 26.27000
C  12.51400 54.21500 25.67000
C  12.80900 54.23100 24.27200
C  14.13200 53.63200 24.11700
C  11.36800 54.71200 26.46900
C  11.91000 54.75300 23.18100
O  12.26900 54.77100 21.97100
C  10.56700 55.35400 23.49600
N  16.86500 52.27500 23.47100
C  15.99700 52.81000 22.61200
C  16.44300 52.50800 21.11400
C  17.81800 51.65600 21.31000
C  17.89300 51.61900 22.83700
C  15.27500 51.72000 20.39100
C  19.11700 52.36100 20.72800
C  19.40800 52.10700 19.22100
N  18.08500 51.15100 25.76500
C  19.11600 50.75400 24.87300
C  20.19800 50.00100 25.64600
C  19.78200 49.96200 26.97300
C  18.52400 50.75100 27.00100
C  21.48500 49.44800 25.11900
C  20.10700 49.55500 28.32000
O  21.06000 48.92400 28.78800
C  19.05800 50.30100 29.28500
C  18.45300 49.37700 30.31700
O  17.48000 48.69900 30.12500
O  19.02600 49.59100 31.59100
C  18.26600 49.12100 32.75900
H  13.19300 53.70900 28.41000
H  14.23600 53.78000 22.01200
H  19.73200 50.45400 22.91400
H  16.24800 51.01000 30.66300
H  14.63800 53.45400 30.28900
H  13.10300 52.02100 30.82500
H  14.20000 50.58000 30.43900
H  13.06200 51.39100 29.27500
H  17.73700 53.44800 30.00000
H  18.21600 52.16100 31.03200
H  16.46200 54.39100 31.85700
H  17.72400 53.56200 32.81800
H  10.55700 54.06200 26.14100
H  11.17100 55.77400 26.32500
H  11.38200 54.60100 27.55400
H  10.71400 56.25300 24.09300
H  10.12800 54.49800 24.01000
H  9.92000 55.56900 22.64600
H  16.62800 53.52500 20.76600
H  17.77900 50.58400 21.11200
H  14.41400 51.50600 21.02400
H  15.64500 50.76000 20.03000
H  14.90700 52.39400 19.61800
H  19.98300 52.23400 21.37700
H  18.92400 53.43100 20.79800
H  19.48300 53.08600 18.74900
H  18.62700 51.45900 18.82300
H  20.36100 51.60000 19.06800
H  22.27100 50.10500 25.49200
H  21.46400 49.41600 24.02900
H  21.78400 48.45300 25.44800
H  19.61000 51.11500 29.75400
H  18.44700 48.04700 32.79300
H  17.22900 49.45300 32.70500
H  18.68800 49.63700 33.62100


