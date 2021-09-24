%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 35.26200 1.75200 29.96200
C  32.85500 2.58300 32.32000
C  37.64200 1.80100 32.35100
C  37.51600 1.12500 27.56800
C  32.74800 1.88400 27.52800
N  35.16600 2.13300 32.11400
C  34.11900 2.35200 32.87300
C  34.63200 2.45300 34.38300
C  36.03300 2.95200 34.09900
C  36.35200 2.21400 32.76000
C  36.05200 4.49500 34.04900
C  34.37500 1.08800 35.15500
C  34.37600 1.20900 36.68200
H  34.84400 2.50900 37.31300
N  37.24000 1.32700 29.98300
C  38.04600 1.32500 31.11900
C  39.36800 0.74800 30.73100
C  39.37200 0.59900 29.27900
C  38.02100 1.10700 28.92100
C  40.46600 0.47300 31.73100
C  40.47100 0.00700 28.28100
O  40.29400 -0.12300 27.06500
C  41.78300 -0.42500 28.90200
N  35.17700 1.65000 27.89700
C  36.22800 1.39600 27.12200
C  35.89000 1.47100 25.60100
C  34.35600 1.46600 25.66700
C  34.03700 1.61900 27.12100
C  36.49900 2.69400 24.86000
C  33.56700 0.32100 24.88900
C  33.00100 -0.79800 25.70000
N  33.26100 2.34500 29.84600
C  32.36300 2.28700 28.83400
C  31.06300 2.71100 29.21500
C  31.19600 2.84400 30.63300
C  32.52400 2.60100 30.90000
C  29.78600 2.98300 28.46300
C  30.54600 3.13100 31.87700
O  29.40100 3.42500 32.14500
C  31.52900 2.83000 33.05300
C  31.55800 3.97100 33.97300
O  31.92600 5.08300 33.62600
O  30.95500 3.71000 35.23000
C  30.85200 4.85600 36.14700
H  38.39100 1.89200 33.13900
H  38.28200 0.96200 26.80800
H  32.00200 1.96400 26.73500
H  34.02600 3.17400 34.93000
H  36.78800 2.59500 34.79900
H  36.64600 4.74400 34.92900
H  35.05700 4.93600 34.12200
H  36.63900 4.75700 33.16900
H  35.14200 0.35900 34.89200
H  33.36400 0.79400 34.87200
H  34.98900 0.39000 37.05700
H  33.35800 1.00500 37.01300
H  40.05900 0.17200 32.69700
H  41.18400 1.27200 31.91300
H  41.03100 -0.39400 31.38900
H  42.08900 0.44800 29.47800
H  42.58200 -0.61900 28.18700
H  41.63400 -1.33400 29.48500
H  36.25600 0.59300 25.06900
H  34.05000 2.44500 25.30000
H  35.69900 3.06000 24.21600
H  37.40400 2.52000 24.28000
H  36.80800 3.41800 25.61400
H  34.19600 -0.03500 24.07300
H  32.66800 0.75100 24.44800
H  31.91500 -0.71000 25.66100
H  33.14600 -0.73200 26.77900
H  33.24300 -1.79400 25.32900
H  29.91400 2.90200 27.38300
H  29.37400 3.93500 28.79900
H  28.98400 2.25500 28.58700
H  31.23100 1.94400 33.61300
H  30.55000 4.60700 37.16400
H  30.10200 5.48600 35.66800
H  31.81400 5.36300 36.08200
Mg 25.47400 50.80800 26.36000
C  23.42800 51.77500 29.14000
C  27.78000 49.80000 28.46900
C  27.48600 50.74800 23.72700
C  22.89100 52.35600 24.24800
N  25.55500 50.75400 28.55700
C  24.64000 51.22200 29.52100
C  25.12200 51.01800 31.00700
C  26.56400 50.28500 30.68300
C  26.67300 50.24400 29.16200
C  27.79400 51.06900 31.31300
C  24.14500 50.22000 31.93900
C  24.75600 49.66700 33.22100
H  24.21600 50.30100 34.51600
N  27.35300 50.24000 26.12300
C  28.22400 49.80700 27.09500
C  29.50900 49.42000 26.55200
C  29.38500 49.70200 25.13400
C  27.97500 50.19600 24.94200
C  30.61900 48.81900 27.27800
C  30.38700 49.50900 24.00900
O  30.13300 49.66800 22.84800
C  31.78500 49.18000 24.46100
N  25.24600 51.52100 24.28400
C  26.25500 51.32600 23.39800
C  25.83300 51.76400 21.97000
C  24.32300 52.13900 22.16700
C  24.10500 51.95100 23.65600
C  26.64800 52.90800 21.27600
C  23.28700 51.29000 21.34900
C  23.20200 51.69600 19.86100
N  23.59000 51.86200 26.56100
C  22.66400 52.33500 25.66000
C  21.52300 52.77100 26.36200
C  21.75600 52.64400 27.68000
C  23.06700 52.09200 27.80400
C  20.26700 53.21900 25.62600
C  21.30400 52.86100 29.05800
O  20.28500 53.34300 29.48100
C  22.42600 52.30600 30.05900
C  23.02400 53.48800 30.75300
O  23.93000 54.23000 30.44700
O  22.33300 53.57900 31.94500
C  22.72800 54.66700 32.87100
H  28.62400 49.48900 29.08800
H  28.23800 50.84300 22.94000
H  22.18300 52.67600 23.48000
H  25.36700 52.00000 31.41200
H  26.54900 49.29900 31.14800
H  28.47500 50.35700 31.77900
H  27.60200 51.88000 32.01500
H  28.41300 51.49500 30.52400
H  23.78800 49.47000 31.23200
H  23.21100 50.69500 32.24000
H  25.83400 49.82400 33.19700
H  24.45000 48.62700 33.33100
H  30.63000 47.83900 26.80100
H  30.42700 48.80700 28.35100
H  31.61000 49.25400 27.15000
H  31.91600 48.09800 24.46700
H  32.15000 49.82600 25.26000
H  32.40300 49.70200 23.73100
H  25.75100 50.83500 21.40500
H  24.20500 53.21400 22.02900
H  26.08400 53.79400 20.98300
H  26.95200 52.49600 20.31400
H  27.50400 53.28000 21.83900
H  22.29300 51.43100 21.77200
H  23.59200 50.25400 21.50300
H  23.57200 52.71200 19.72200
H  22.15600 51.57400 19.58000
H  23.82700 51.09700 19.19800
H  19.37700 53.06800 26.23700
H  20.20500 52.68100 24.68000
H  20.36300 54.25800 25.31000
H  22.01700 51.57800 30.76000
H  21.79500 55.04700 33.28700
H  23.27800 55.48100 32.39900
H  23.31400 54.08500 33.58300


