%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(current)

TDDFT excited states in gaussian

0 1
Mg 35.66000 1.74600 29.23200
C  33.47200 2.34100 31.89500
C  38.12100 1.24100 31.62500
C  37.79000 1.61500 26.77900
C  33.01700 2.53400 27.01000
N  35.87600 2.18500 31.54100
C  34.82500 2.24700 32.36500
C  35.34000 2.18200 33.76400
C  36.88000 2.27300 33.63100
C  37.00600 1.84600 32.18200
C  37.59500 3.59000 33.93800
C  34.89600 0.87200 34.48900
C  34.26200 0.97200 35.91400
H  34.68800 2.14100 36.85600
N  37.69800 1.31100 29.19900
C  38.45500 0.95400 30.30700
C  39.74800 0.60100 29.81800
C  39.73600 0.82500 28.38300
C  38.39000 1.25600 28.03600
C  40.94400 0.20000 30.70900
C  40.89400 0.72300 27.34900
O  40.78100 0.94900 26.18600
C  42.22500 0.23800 27.88900
N  35.39600 2.02600 27.20700
C  36.47200 1.98000 26.39000
C  36.13800 2.22700 24.89300
C  34.55100 2.25100 24.96900
C  34.26500 2.24500 26.46600
C  36.85300 3.59900 24.36100
C  33.86400 1.14100 24.11300
C  33.41700 -0.12600 24.82000
N  33.60700 2.29400 29.34300
C  32.66300 2.48800 28.37100
C  31.34700 2.76700 29.04600
C  31.64200 2.71200 30.42700
C  33.04500 2.46100 30.58700
C  30.05700 3.20300 28.40600
C  31.07900 2.83900 31.76900
O  29.91300 3.09600 32.08200
C  32.20100 2.44500 32.72500
C  32.26400 3.54000 33.75100
O  32.77700 4.63100 33.58200
O  31.49900 3.13900 34.75300
C  31.03100 4.13600 35.70300
H  38.92600 1.22300 32.36200
H  38.44500 1.39800 25.93300
H  32.19900 2.76200 26.32300
H  34.93300 3.02400 34.32300
H  37.28500 1.57100 34.36000
H  36.96900 4.41100 34.28800
H  38.13600 3.93900 33.05800
H  38.32400 3.36100 34.71500
H  35.66700 0.10200 34.50700
H  34.14500 0.53000 33.77700
H  34.54900 0.03000 36.38300
H  33.17800 1.00900 35.80900
H  41.84000 0.79100 30.51900
H  41.29800 -0.80900 30.49700
H  40.73900 0.22800 31.78000
H  42.13800 -0.68600 28.46100
H  42.66900 1.05000 28.46500
H  42.94900 0.14100 27.08000
H  36.48000 1.41700 24.24900
H  34.11700 3.15500 24.54200
H  36.02900 4.30300 24.24400
H  37.25500 3.51600 23.35100
H  37.53000 4.05700 25.08300
H  34.65800 0.79900 23.44900
H  33.05000 1.52900 23.50100
H  34.25200 -0.36700 25.47900
H  33.30700 -0.80900 23.97800
H  32.56000 -0.11800 25.49300
H  30.04900 3.41600 27.33700
H  29.71100 4.11800 28.88700
H  29.33300 2.43700 28.68200
H  31.98900 1.57300 33.34400
H  31.22900 5.19300 35.52800
H  31.60100 3.89200 36.60000
H  30.07700 3.79700 36.10700


