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
Mg 40.74100 8.40300 29.69900
C  42.66000 10.02500 32.01600
C  38.73200 7.21500 32.07700
C  39.17900 6.75200 27.36300
C  43.28300 9.22800 27.34100
N  40.79300 8.45600 31.86800
C  41.50500 9.39200 32.54400
C  40.94800 9.63500 33.92000
C  39.85300 8.49700 34.00200
C  39.82700 7.94800 32.59200
C  40.10000 7.44100 35.10600
C  40.44100 11.18100 34.18500
C  40.87000 11.88200 35.48300
H  41.23900 10.96600 36.69200
N  39.13800 7.34400 29.69200
C  38.39900 6.91700 30.74900
C  37.34700 6.00300 30.37600
C  37.33700 6.02300 29.01400
C  38.65500 6.62800 28.64700
C  36.61000 5.29800 31.50000
C  36.33400 5.36600 28.04500
O  36.44200 5.38700 26.84100
C  35.20700 4.65100 28.73000
N  41.30000 7.91900 27.67200
C  40.36700 7.31400 26.91500
C  40.80200 7.26400 25.44700
C  41.95100 8.32700 25.38700
C  42.23200 8.49800 26.89400
C  41.21600 5.86100 24.82900
C  41.58000 9.72400 24.76100
C  40.48400 10.59900 25.31000
N  42.57200 9.48100 29.57000
C  43.55300 9.66900 28.59100
C  44.67900 10.46300 29.09900
C  44.36400 10.63700 30.45500
C  43.08700 10.02200 30.64800
C  45.85600 11.01300 28.33500
C  44.78800 11.17000 31.66800
O  45.80500 11.75000 32.01500
C  43.61600 11.01200 32.67100
C  44.29700 10.49000 33.87800
O  44.63700 9.31800 34.05700
O  44.61300 11.49900 34.74100
C  45.54500 11.24000 35.84300
H  37.99000 7.00700 32.85100
H  38.68100 6.27000 26.52000
H  43.95700 9.44200 26.50800
H  41.68000 9.37000 34.68300
H  38.87700 8.95900 34.14900
H  40.48500 6.53200 34.64400
H  39.14600 7.20600 35.57800
H  40.83000 7.87400 35.79000
H  39.35300 11.12200 34.15300
H  40.80700 11.84200 33.40000
H  40.04300 12.46600 35.88600
H  41.65400 12.60800 35.26800
H  37.22900 4.66200 32.13300
H  35.89000 4.51000 31.28000
H  36.12800 6.14400 31.99000
H  34.67500 5.32500 29.40100
H  35.72300 3.80100 29.17600
H  34.37800 4.27900 28.12900
H  39.88500 7.56300 24.93800
H  42.83200 7.86300 24.94300
H  41.11200 5.08600 25.58800
H  42.28100 5.83900 24.60100
H  40.70000 5.56000 23.91600
H  41.36300 9.40900 23.74000
H  42.43100 10.40400 24.73900
H  39.92000 10.04700 26.06300
H  39.77700 10.86400 24.52400
H  40.79900 11.50000 25.83600
H  46.78200 10.48700 28.56300
H  45.97200 12.07500 28.55200
H  45.60900 10.77100 27.30100
H  43.07400 11.94700 32.81400
H  46.05900 12.14500 36.16900
H  46.18500 10.38200 35.64300
H  45.04400 10.86200 36.73400


