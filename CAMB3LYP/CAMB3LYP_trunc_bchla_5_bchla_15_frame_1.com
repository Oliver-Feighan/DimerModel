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
Mg 46.36700 34.87300 28.25800
C  44.78600 32.87200 30.57800
C  46.38600 37.42000 30.77600
C  47.55600 36.98300 26.07000
C  45.97900 32.44300 25.82700
N  45.44200 35.14400 30.41100
C  45.07300 34.11700 31.15800
C  45.12500 34.48000 32.66900
C  45.08600 36.04500 32.56800
C  45.67500 36.28100 31.16800
C  43.63100 36.61800 32.70900
C  46.46000 33.91200 33.20000
C  46.69200 34.00900 34.71400
H  45.77100 34.85100 35.57400
N  47.01100 36.84800 28.33900
C  46.96000 37.70200 29.49400
C  47.69800 38.90500 29.13400
C  48.14000 38.73300 27.78500
C  47.54200 37.48800 27.30500
C  47.97900 40.11100 30.10600
C  49.09900 39.67100 26.98000
O  49.27000 39.43700 25.81500
C  49.74700 40.84400 27.64100
N  46.80600 34.69600 26.23000
C  47.30800 35.72100 25.50500
C  47.46400 35.39100 24.02100
C  47.08600 33.89200 23.98700
C  46.58100 33.65400 25.46800
C  46.60500 36.26500 23.03500
C  48.11600 32.92900 23.45400
C  49.59100 32.99000 23.91900
N  45.55600 33.03000 28.12500
C  45.46900 32.14000 27.09300
C  44.91700 30.91500 27.63200
C  44.64100 31.10300 29.03000
C  45.05100 32.44500 29.25900
C  44.67700 29.67800 26.82400
C  44.09100 30.54400 30.21900
O  43.56800 29.43800 30.47700
C  43.93800 31.83400 31.19500
C  44.25300 31.34600 32.61000
O  45.22200 30.77300 32.93400
O  43.18900 31.65700 33.40800
C  43.38400 31.62900 34.84300
H  46.32800 38.21100 31.52600
H  47.90200 37.66100 25.28700
H  45.81100 31.64300 25.10400
H  44.23500 34.14900 33.20400
H  45.76700 36.36500 33.35700
H  43.51000 37.05600 33.69900
H  42.93300 35.79700 32.54400
H  43.41600 37.43400 32.01900
H  47.26100 34.50400 32.75800
H  46.66000 32.91100 32.81900
H  47.74400 34.11200 34.98100
H  46.57400 32.98400 35.06600
H  47.70900 41.07100 29.66500
H  49.03500 40.03600 30.36500
H  47.45400 39.94800 31.04800
H  50.66300 40.94800 27.05900
H  50.08200 40.62400 28.65400
H  49.16800 41.73600 27.88000
H  48.52200 35.54500 23.81000
H  46.16900 33.69700 23.43100
H  46.31700 35.59300 22.22700
H  47.32700 37.01000 22.70000
H  45.80600 36.77400 23.57400
H  48.21500 33.30700 22.43700
H  47.86400 31.86800 23.45400
H  50.35800 33.20100 23.17500
H  49.94100 32.01500 24.25800
H  49.69200 33.75700 24.68700
H  44.20400 28.97200 27.50600
H  45.62900 29.23800 26.52400
H  44.15300 29.83100 25.88100
H  42.91200 32.17200 31.04800
H  42.52900 30.96500 34.97100
H  43.38800 32.61600 35.30600
H  44.28100 31.10000 35.16600


