%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

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
Mg -2.59000 34.01700 26.79400
C  -3.62000 32.41400 29.77300
C  -1.15300 36.45800 28.71900
C  -2.20600 35.82500 24.04600
C  -4.51700 31.70900 25.00200
N  -2.55200 34.47200 29.00300
C  -2.95200 33.59300 30.01300
C  -2.37700 34.05700 31.37400
C  -1.60200 35.37500 31.00700
C  -1.71600 35.46600 29.46000
C  -1.91400 36.61400 31.82800
C  -1.50900 33.04200 32.13400
C  -1.98900 32.76400 33.58500
H  -0.95800 32.99700 34.66700
N  -1.70200 35.86000 26.44900
C  -1.19200 36.74900 27.34900
C  -0.55100 37.80700 26.60500
C  -0.98500 37.71000 25.21400
C  -1.68000 36.43200 25.23200
C  0.22600 39.02700 27.32500
C  -1.03600 38.70200 24.13300
O  -1.40600 38.53400 23.00600
C  -0.57600 40.11000 24.49900
N  -3.46600 33.87800 24.82300
C  -3.08300 34.73100 23.85300
C  -3.52500 34.22200 22.47900
C  -3.92700 32.76700 22.76900
C  -3.94500 32.73700 24.28100
C  -4.74000 35.08300 21.90100
C  -2.77000 31.79500 22.20500
C  -1.44300 31.78800 23.01200
N  -3.75200 32.29300 27.19100
C  -4.44100 31.48200 26.37200
C  -5.05500 30.45500 27.17000
C  -4.72100 30.75200 28.47800
C  -3.96900 31.93000 28.47700
C  -5.92100 29.29000 26.73000
C  -4.97000 30.27500 29.84100
O  -5.63000 29.44100 30.39400
C  -4.20300 31.28400 30.73600
C  -5.22200 31.86200 31.64700
O  -6.22400 32.38000 31.22200
O  -4.94100 31.66200 32.98900
C  -6.05600 31.94700 33.95300
H  -0.52600 37.10900 29.33200
H  -2.00600 36.32000 23.09300
H  -4.95700 30.83400 24.51900
H  -3.20300 34.37600 32.01000
H  -0.55500 35.23200 31.27500
H  -2.84600 36.49700 32.37900
H  -1.95200 37.53200 31.24100
H  -1.13500 36.66900 32.58800
H  -0.45800 33.28400 32.29200
H  -1.47400 32.09100 31.60300
H  -2.46500 31.78700 33.66800
H  -2.77300 33.48500 33.81300
H  1.18600 39.15100 26.82400
H  0.47900 38.88500 28.37500
H  -0.26700 39.99700 27.38900
H  0.46900 39.96700 24.77500
H  -1.19700 40.39100 25.34900
H  -0.59300 40.82400 23.67500
H  -2.71600 34.17900 21.74900
H  -4.90200 32.39200 22.45700
H  -4.66500 35.04300 20.81400
H  -4.79700 36.14200 22.15500
H  -5.62600 34.53900 22.22800
H  -2.53100 31.98600 21.15900
H  -3.07000 30.75600 22.34200
H  -0.61000 31.26900 22.53900
H  -1.63900 31.28300 23.95800
H  -1.17900 32.82600 23.21700
H  -5.41600 28.57400 26.08200
H  -6.72000 29.76700 26.16300
H  -6.40800 28.78800 27.56600
H  -3.52700 30.64600 31.30500
H  -5.48500 32.48600 34.70900
H  -6.48400 31.04500 34.39000
H  -6.92400 32.52400 33.63500


