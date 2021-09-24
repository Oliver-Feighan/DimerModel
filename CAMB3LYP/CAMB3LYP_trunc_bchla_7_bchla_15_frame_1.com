%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 26.07000 0.32300 29.17500
C  27.93000 0.00300 32.17200
C  23.23700 0.37300 31.22200
C  24.20000 0.73500 26.50800
C  28.89900 -0.55900 27.35600
N  25.68800 0.23300 31.42600
C  26.58700 0.25300 32.44900
C  25.99300 0.25100 33.86800
C  24.43500 0.24800 33.45400
C  24.48700 0.27000 31.95500
C  23.63400 -1.02300 33.94300
C  26.39900 1.50900 34.70300
C  27.30600 1.47800 35.92700
H  27.07500 2.41400 37.10600
N  24.00400 0.60600 28.87400
C  22.98600 0.52400 29.83500
C  21.72400 0.61300 29.13500
C  21.97300 0.88000 27.79100
C  23.44300 0.80400 27.68300
C  20.46800 0.64800 29.93700
C  21.02800 1.08200 26.60700
O  21.39200 1.28100 25.44400
C  19.58500 1.18400 26.91000
N  26.42900 0.02500 27.20000
C  25.57000 0.46200 26.23000
C  26.30100 0.60400 24.85200
C  27.70600 0.10100 25.15700
C  27.68000 -0.24400 26.67200
C  25.55900 -0.10800 23.67400
C  28.73300 1.18300 24.79000
C  29.96300 0.77300 24.06500
N  28.00500 -0.28200 29.61000
C  29.05600 -0.58000 28.78000
C  30.24500 -0.86800 29.52500
C  29.92400 -0.56000 30.93700
C  28.52000 -0.19800 30.89800
C  31.58600 -1.33700 28.93500
C  30.40700 -0.48700 32.33900
O  31.55800 -0.58600 32.79700
C  29.11500 -0.11100 33.17500
C  29.05800 -1.28400 34.10300
O  28.53100 -2.37500 33.92200
O  29.56000 -0.97800 35.26800
C  29.75400 -2.03800 36.28700
H  22.46700 0.52900 31.98000
H  23.69600 0.81600 25.54300
H  29.79300 -0.76400 26.76400
H  26.26300 -0.69600 34.33700
H  23.82700 1.11200 33.72400
H  24.21800 -1.88400 34.26700
H  22.99700 -1.39000 33.13900
H  22.99100 -0.72400 34.77100
H  25.46900 1.94600 35.06700
H  26.90100 2.23700 34.06700
H  28.23400 1.84700 35.49000
H  27.27900 0.42100 36.19200
H  19.92700 1.58400 29.80500
H  20.67800 0.53600 31.00100
H  19.76300 -0.13800 29.66800
H  19.43300 1.90200 27.71600
H  19.18900 0.27500 27.36300
H  18.97400 1.53500 26.07900
H  26.25300 1.68100 24.69200
H  27.96600 -0.74800 24.52500
H  24.56100 -0.46700 23.92500
H  26.08600 -0.97800 23.28400
H  25.56500 0.53500 22.79400
H  29.00400 1.58600 25.76500
H  28.34900 1.97700 24.15000
H  30.78000 1.14300 24.68400
H  29.93300 1.13400 23.03700
H  30.14200 -0.30200 24.03500
H  31.60000 -2.42600 28.91000
H  32.38900 -0.93100 29.55000
H  31.71700 -1.01400 27.90200
H  29.35100 0.80400 33.71800
H  28.83300 -2.53300 36.59300
H  30.27100 -1.70300 37.18600
H  30.41000 -2.86500 36.01700
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


