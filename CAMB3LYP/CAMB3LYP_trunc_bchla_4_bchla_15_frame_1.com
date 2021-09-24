%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 8.82800 2.97900 28.39100
C  10.36500 1.60400 31.10700
C  7.59700 5.54500 30.49200
C  7.11400 4.23300 25.80200
C  10.13400 0.47400 26.43700
N  8.69300 3.35100 30.57800
C  9.52500 2.68500 31.48700
C  9.40600 3.27200 32.92200
C  8.46000 4.46400 32.61800
C  8.16100 4.43600 31.13000
C  7.12300 4.69600 33.44300
C  10.76000 3.72700 33.62500
C  10.50900 4.32800 35.03400
H  10.96800 3.41400 36.17100
N  7.54400 4.72800 28.17300
C  7.38800 5.73600 29.10000
C  6.47100 6.69300 28.49200
C  6.33200 6.29300 27.12400
C  7.02200 5.05300 26.95500
C  5.95900 7.85000 29.26000
C  5.41800 7.05100 26.09900
O  5.12000 6.60100 24.97400
C  4.78500 8.36300 26.48300
N  8.63500 2.44800 26.51900
C  7.83100 3.12600 25.56000
C  7.67200 2.29800 24.30500
C  8.94400 1.45600 24.38400
C  9.32300 1.47600 25.85800
C  6.36700 1.43300 24.30400
C  10.11600 2.06600 23.58300
C  10.22200 1.75100 22.11500
N  9.99100 1.29200 28.61900
C  10.57300 0.45400 27.74600
C  11.35400 -0.51100 28.46200
C  11.33200 -0.10700 29.78000
C  10.51600 1.03400 29.82700
C  12.08500 -1.61500 27.78100
C  11.82500 -0.42000 31.05300
O  12.58200 -1.29100 31.43600
C  11.20500 0.70900 32.00300
C  10.44000 0.02800 33.08200
O  9.42100 -0.60800 32.96400
O  10.93800 0.32100 34.28700
C  10.21200 -0.27100 35.47600
H  7.38800 6.37900 31.16500
H  6.56600 4.63400 24.94700
H  10.68600 -0.27500 25.86600
H  8.89500 2.54300 33.55100
H  9.00500 5.40000 32.74400
H  6.85500 5.73300 33.64700
H  7.10100 4.02900 34.30500
H  6.33500 4.37900 32.76000
H  11.20100 4.41000 32.90000
H  11.41200 2.85700 33.70200
H  9.47300 4.65600 35.11900
H  11.13500 5.21700 35.10500
H  4.89200 7.69700 29.42300
H  6.27900 8.73500 28.71000
H  6.46700 7.93200 30.22100
H  4.23800 8.77300 25.63400
H  5.63600 9.03200 26.61800
H  4.16200 8.37000 27.37700
H  7.64200 2.93000 23.41700
H  8.80000 0.42300 24.06800
H  6.45900 0.58300 24.97900
H  6.13800 1.09500 23.29300
H  5.57600 2.10700 24.63100
H  11.08800 1.85700 24.02800
H  10.10000 3.15100 23.68700
H  10.17300 2.70400 21.58900
H  9.27100 1.26300 21.90100
H  11.07000 1.12100 21.84600
H  13.12700 -1.33900 27.62400
H  11.63900 -2.06200 26.89200
H  12.10800 -2.42700 28.50800
H  12.02700 1.28100 32.43400
H  9.15600 -0.02800 35.35600
H  10.52300 0.16000 36.42800
H  10.37400 -1.34900 35.46000
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


