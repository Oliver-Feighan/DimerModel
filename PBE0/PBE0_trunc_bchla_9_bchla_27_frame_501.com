%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

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
Mg -5.09800 24.63300 26.68400
C  -3.70200 26.36700 29.52200
C  -6.11300 22.14300 28.93500
C  -6.48000 23.09600 24.13300
C  -3.72900 27.06700 24.58100
N  -4.87100 24.33000 29.01500
C  -4.47600 25.23800 29.94600
C  -4.87000 24.81900 31.35300
C  -5.34300 23.26900 31.18300
C  -5.42800 23.22100 29.63000
C  -4.42000 22.22500 31.86000
C  -5.97200 25.74300 31.97300
C  -5.56400 26.39400 33.32000
H  -6.48200 25.91600 34.58700
N  -6.12400 22.92100 26.55900
C  -6.45400 21.95500 27.56900
C  -7.23800 20.80800 27.05300
C  -7.37200 21.15200 25.66800
C  -6.65700 22.41400 25.41400
C  -7.70400 19.61500 27.89000
C  -8.01000 20.37100 24.58300
O  -8.30600 20.80000 23.47300
C  -8.47200 18.93600 24.90400
N  -5.23600 25.12700 24.71500
C  -5.84300 24.30800 23.82600
C  -5.57400 24.69000 22.39700
C  -4.73400 25.98700 22.49800
C  -4.58300 26.13300 24.02900
C  -4.97900 23.58700 21.54700
C  -5.55400 27.25800 21.86300
C  -4.76800 28.10300 20.72000
N  -3.71900 26.23300 26.92800
C  -3.33300 27.18000 25.94900
C  -2.55400 28.22500 26.62300
C  -2.81300 27.99000 27.99900
C  -3.40400 26.73200 28.12600
C  -1.71900 29.35800 26.03400
C  -2.39200 28.46300 29.34800
O  -1.66000 29.39300 29.78700
C  -3.08100 27.46800 30.36000
C  -2.03500 26.94900 31.30200
O  -0.94600 26.57900 30.96500
O  -2.47500 27.07900 32.59000
C  -1.38400 26.91100 33.58300
H  -6.28600 21.33800 29.65200
H  -6.91400 22.63300 23.24400
H  -3.32400 27.85900 23.94800
H  -3.94900 24.82500 31.93600
H  -6.36100 23.10200 31.53600
H  -3.94100 21.58900 31.11500
H  -4.98800 21.61200 32.56100
H  -3.64700 22.62700 32.51500
H  -6.89500 25.19600 32.16700
H  -6.22100 26.57200 31.31000
H  -5.77000 27.46400 33.35800
H  -4.50500 26.32600 33.56700
H  -7.71300 19.74200 28.97300
H  -7.20000 18.67100 27.68300
H  -8.74200 19.45800 27.59500
H  -7.62000 18.31700 25.18500
H  -9.03200 18.48700 24.08300
H  -9.15600 18.98900 25.75100
H  -6.53500 24.86100 21.91000
H  -3.75100 25.85800 22.04600
H  -5.68000 23.16900 20.82400
H  -4.66700 22.69800 22.09500
H  -4.09700 23.98900 21.04800
H  -5.86900 27.99500 22.60100
H  -6.48500 26.95600 21.38400
H  -4.00100 27.42900 20.33800
H  -4.32800 29.01700 21.12100
H  -5.40200 28.28200 19.85200
H  -1.71900 30.24300 26.67000
H  -2.23000 29.70200 25.13400
H  -0.77900 28.87000 25.77600
H  -3.74400 28.08300 30.96800
H  -0.90300 25.95500 33.37700
H  -1.91800 26.79500 34.52600
H  -0.58100 27.64700 33.62200


