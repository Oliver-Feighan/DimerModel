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
Mg 29.46400 59.14500 41.03800
C  26.57800 57.38900 40.28900
C  31.34900 56.53100 39.78200
C  32.15000 61.18900 41.08600
C  27.35400 61.82000 41.94900
N  28.99200 57.29200 39.87100
C  27.77200 56.74800 39.89300
C  27.77900 55.38900 39.12700
C  29.32900 54.99700 39.21800
C  29.96700 56.34400 39.73200
C  29.58100 53.72700 40.06600
C  27.37400 55.59400 37.68500
C  26.46000 54.52700 37.10300
H  26.57800 54.02100 35.61600
N  31.53300 58.78800 40.65000
C  32.11000 57.71400 40.08400
C  33.51500 57.95800 40.03300
C  33.80600 59.24700 40.57400
C  32.47900 59.81600 40.84000
C  34.47700 56.89100 39.66100
C  35.16200 59.85800 40.59300
O  36.13000 59.21400 40.22400
C  35.39700 61.26400 41.13200
N  29.71300 61.28200 41.56500
C  30.94400 61.84600 41.41400
C  30.85900 63.33400 41.52000
C  29.29200 63.52200 41.79300
C  28.69900 62.13600 41.84900
C  31.79900 63.96200 42.56600
C  28.51800 64.48000 40.78300
C  27.64200 65.63900 41.31400
N  27.43700 59.54700 41.14000
C  26.74400 60.62700 41.58900
C  25.41200 60.28000 41.81900
C  25.25400 58.98400 41.34800
C  26.48900 58.61100 40.88000
C  24.32700 61.25500 42.36300
C  24.37200 57.92800 41.16800
O  23.15900 57.91800 41.34600
C  25.14900 56.75500 40.39700
C  25.25800 55.58200 41.18700
O  25.99300 55.31600 42.14100
O  24.33400 54.68200 40.66500
C  24.17500 53.35900 41.24400
H  31.95200 55.65500 39.53500
H  32.91900 61.96200 41.15300
H  26.80500 62.68200 42.33400
H  27.08300 54.71900 39.63100
H  29.81700 54.88800 38.25000
H  28.66900 53.37400 40.54700
H  30.34000 53.93800 40.82000
H  30.00600 52.97200 39.40500
H  28.16800 55.74800 36.95400
H  26.92700 56.58800 37.67000
H  25.44000 54.85100 37.30900
H  26.66400 53.69100 37.77200
H  34.81300 56.52500 40.63100
H  35.29700 57.28300 39.05900
H  34.02200 56.08300 39.08800
H  35.10600 61.93900 40.32700
H  36.47900 61.37300 41.20200
H  34.84400 61.42100 42.05800
H  31.14300 63.61000 40.50500
H  29.10900 64.04800 42.73000
H  32.59600 64.44300 41.99800
H  32.19100 63.23400 43.27700
H  31.26900 64.75100 43.09900
H  27.85100 63.86700 40.17700
H  29.18800 65.01800 40.11200
H  26.69000 65.11600 41.40400
H  27.54400 66.41700 40.55800
H  27.77900 66.00900 42.33000
H  24.06200 61.95200 41.56900
H  24.63400 61.90700 43.18100
H  23.42900 60.71200 42.65900
H  24.61600 56.61800 39.45600
H  23.13300 53.07400 41.10200
H  24.52900 53.25500 42.27000
H  24.72400 52.65900 40.61400


