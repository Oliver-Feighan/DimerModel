%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 24.45900 -7.80400 45.73000
C  26.58600 -5.04500 44.57700
C  21.88700 -6.07400 44.19000
C  22.84100 -10.55200 45.75100
C  27.52200 -9.37000 46.76800
N  24.34300 -5.76600 44.35600
C  25.30500 -4.83000 44.04500
C  24.68800 -3.74900 43.08400
C  23.16700 -3.94300 43.21100
C  23.11900 -5.33100 43.95600
C  22.33200 -2.87300 44.03000
C  25.23700 -3.82200 41.61800
C  25.48800 -2.49900 40.92700
H  24.41700 -2.03400 39.89100
N  22.60400 -8.23500 45.19100
C  21.64700 -7.40100 44.68000
C  20.38100 -8.15500 44.46900
C  20.67900 -9.47200 44.90600
C  22.09500 -9.47000 45.26500
C  19.08200 -7.63000 43.83600
C  19.81800 -10.73300 44.78900
O  18.67600 -10.59500 44.35900
C  20.22700 -12.02800 45.16700
N  25.20700 -9.77600 45.95500
C  24.19400 -10.69000 46.16700
C  24.84300 -12.02000 46.71400
C  26.38100 -11.70100 46.72900
C  26.35300 -10.16200 46.56900
C  24.30300 -12.52000 48.08200
C  27.32200 -12.42400 45.68300
C  27.76100 -13.86800 45.98200
N  26.63200 -7.37200 45.76700
C  27.69000 -8.05400 46.33400
C  28.89800 -7.26500 46.33900
C  28.50100 -6.05100 45.66300
C  27.12300 -6.13300 45.25800
C  30.21600 -7.74700 46.86200
C  28.99300 -4.84100 45.17500
O  30.13700 -4.41700 45.08200
C  27.84800 -4.18900 44.34700
C  27.67100 -2.81900 44.95900
O  27.62100 -1.78000 44.29600
O  27.39300 -2.80300 46.38500
C  26.88100 -1.51100 46.83900
H  21.04400 -5.51700 43.77700
H  22.37000 -11.50900 45.98500
H  28.34500 -9.94300 47.20000
H  24.83100 -2.76600 43.53500
H  22.69900 -4.12200 42.24300
H  21.79800 -2.30300 43.27000
H  22.98800 -2.17400 44.54800
H  21.67600 -3.27700 44.80100
H  24.55300 -4.27200 40.89800
H  26.19500 -4.34200 41.59800
H  26.49400 -2.54600 40.51000
H  25.50300 -1.72600 41.69600
H  18.22200 -8.02700 44.37600
H  19.15100 -7.95200 42.79700
H  19.01200 -6.54700 43.92800
H  20.35400 -12.03800 46.25000
H  21.09800 -12.40800 44.63300
H  19.40000 -12.62300 44.77800
H  24.71500 -12.87500 46.05100
H  26.80800 -11.85300 47.72000
H  25.17900 -12.58000 48.72800
H  23.87600 -13.51200 47.93800
H  23.57100 -11.76000 48.35800
H  28.20300 -11.78200 45.65400
H  26.82800 -12.52100 44.71600
H  27.32600 -14.00800 46.97200
H  28.83900 -13.92900 46.13800
H  27.42200 -14.64300 45.29500
H  30.43000 -8.65700 46.30300
H  30.17700 -8.08000 47.89900
H  31.04400 -7.08900 46.59800
H  28.22500 -4.06500 43.33200
H  26.34600 -1.90700 47.70300
H  26.35400 -0.97600 46.04900
H  27.73900 -0.96400 47.23000
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


