%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
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
Mg 35.49600 50.25400 25.27000
C  35.17700 48.23300 28.13300
C  33.95100 52.80200 27.06200
C  35.52400 52.00900 22.45100
C  36.25500 47.31600 23.35100
N  34.65000 50.47900 27.32000
C  34.75900 49.53200 28.34800
C  34.55800 50.22300 29.68400
C  33.90100 51.58900 29.33600
C  34.18900 51.63300 27.80500
C  32.41200 51.77300 29.76700
C  35.85100 50.33500 30.56400
C  35.57700 50.30300 32.13200
H  34.12200 50.22500 32.60200
N  34.96300 52.12400 24.87400
C  34.41800 53.04500 25.75200
C  34.29200 54.33000 25.12700
C  34.61900 54.10800 23.75100
C  35.08400 52.73500 23.66700
C  33.72300 55.59000 25.84800
C  34.52600 55.05800 22.60200
O  34.69900 54.67600 21.42000
C  34.14200 56.56700 22.60300
N  35.57000 49.71100 23.15300
C  35.66800 50.65100 22.19400
C  36.16600 50.02300 20.86800
C  36.37900 48.52900 21.13700
C  36.08600 48.52400 22.70200
C  35.30200 50.41100 19.70700
C  37.75400 47.93500 20.76300
C  38.97000 48.58500 21.46800
N  35.63500 48.14500 25.59300
C  36.09300 47.16800 24.77600
C  36.32500 45.99400 25.61600
C  36.01700 46.32200 26.96600
C  35.54300 47.65500 26.90400
C  36.93800 44.67500 25.15300
C  35.96700 45.83900 28.30900
O  36.34100 44.80600 28.85800
C  35.32300 47.11200 29.16700
C  34.10500 46.63100 29.93300
O  33.00700 46.49100 29.43100
O  34.39500 46.48000 31.27000
C  33.50200 45.67800 32.10000
H  33.57900 53.56200 27.75200
H  35.60400 52.66100 21.57900
H  36.61900 46.49900 22.72500
H  33.85700 49.64500 30.28600
H  34.49900 52.36000 29.82100
H  31.92700 52.02300 28.82400
H  32.30200 52.65000 30.40600
H  31.87100 50.93900 30.21400
H  36.34800 51.29400 30.42100
H  36.54700 49.55500 30.25800
H  35.89400 51.26500 32.53500
H  36.14900 49.47100 32.54200
H  33.11100 56.13200 25.12600
H  34.50500 56.27400 26.17600
H  33.12800 55.50300 26.75700
H  34.79000 57.12800 23.27700
H  33.09000 56.73000 22.83700
H  34.45500 57.11000 21.71200
H  37.07600 50.56600 20.61600
H  35.49800 48.04500 20.71400
H  35.74600 51.19400 19.09200
H  34.31200 50.69900 20.06100
H  35.07200 49.59000 19.02700
H  37.83100 48.14500 19.69600
H  37.75500 46.85700 20.92700
H  39.44400 47.77500 22.02200
H  38.57400 49.35600 22.12900
H  39.64500 49.11900 20.79900
H  37.67500 44.28600 25.85500
H  37.52900 44.80600 24.24700
H  36.12400 43.98800 24.92000
H  36.09300 47.35400 29.90000
H  32.48700 45.91500 31.78000
H  33.56400 45.88200 33.16900
H  33.81500 44.64500 31.94800


