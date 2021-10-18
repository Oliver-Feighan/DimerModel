%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 25.26700 0.12000 29.48500
C  27.14100 -0.08700 32.36300
C  22.44200 0.65100 31.48400
C  23.33900 0.24400 26.71700
C  27.97000 -0.76800 27.63500
N  24.82700 -0.00300 31.72000
C  25.83000 0.13900 32.68700
C  25.26100 0.57400 34.05200
C  23.73000 0.68200 33.70700
C  23.61100 0.48500 32.21100
C  22.87000 -0.30200 34.49000
C  25.80700 1.87500 34.74000
C  26.64900 1.68700 36.11600
H  25.98700 2.43400 37.27800
N  23.21300 0.29800 29.14100
C  22.23900 0.64000 30.09300
C  20.97400 1.01400 29.45200
C  21.16000 0.82900 27.97700
C  22.63600 0.50000 27.85200
C  19.74900 1.44500 30.22800
C  20.23300 0.96600 26.86300
O  20.59800 0.66100 25.71800
C  18.79000 1.40700 27.03500
N  25.61300 -0.29200 27.55300
C  24.69700 -0.04800 26.51300
C  25.29200 -0.18000 25.05300
C  26.81800 -0.54700 25.43400
C  26.79800 -0.51500 26.97900
C  24.57500 -1.29000 24.13700
C  27.83300 0.46900 24.80600
C  29.13600 -0.06400 24.29800
N  27.18000 -0.44400 29.94400
C  28.21700 -0.73500 29.05300
C  29.45900 -0.85000 29.73800
C  29.08500 -0.49600 31.08700
C  27.69300 -0.34600 31.12700
C  30.82000 -1.11200 29.24500
C  29.58400 -0.41400 32.42400
O  30.70000 -0.50100 32.94800
C  28.28000 -0.28500 33.33500
C  28.18400 -1.48700 34.20500
O  27.52600 -2.44100 34.06000
O  28.92900 -1.30000 35.36700
C  28.48600 -2.07100 36.51000
H  21.58300 0.90400 32.10900
H  22.84400 0.32400 25.74700
H  28.87000 -0.90500 27.03200
H  25.38700 -0.24700 34.75700
H  23.31500 1.64400 34.00900
H  22.43700 -1.01900 33.79300
H  22.09100 0.33900 34.90100
H  23.35600 -0.79100 35.33400
H  24.95600 2.52500 34.94500
H  26.49900 2.40700 34.08700
H  27.58200 2.20000 35.88200
H  26.76100 0.62900 36.35200
H  19.34600 0.54600 30.69600
H  19.02800 1.84900 29.51900
H  19.97500 2.10600 31.06500
H  18.28300 1.10300 26.11900
H  18.67500 2.47500 27.22300
H  18.18400 0.86400 27.76000
H  25.32400 0.80000 24.57800
H  27.16700 -1.50900 25.05800
H  25.34100 -1.79600 23.54900
H  23.91500 -0.92500 23.35100
H  24.00900 -2.02000 24.71600
H  28.05800 1.32200 25.44600
H  27.45300 0.92700 23.89300
H  29.88300 0.35500 24.97200
H  29.54600 0.28700 23.35100
H  29.13600 -1.15300 24.35500
H  30.88400 -0.76300 28.21400
H  30.89900 -2.19900 29.25700
H  31.53400 -0.73600 29.97800
H  28.35300 0.59900 33.96800
H  28.85600 -3.07500 36.30500
H  27.39700 -2.05300 36.47000
H  28.92200 -1.70000 37.43800
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


