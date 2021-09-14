%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 46.81700 25.04400 28.06300
C  47.43500 27.22400 30.70400
C  46.22200 22.54500 30.24300
C  46.68300 22.91500 25.44100
C  47.91200 27.51000 25.91600
N  46.85500 24.93200 30.27600
C  46.92700 25.99900 31.13500
C  46.69200 25.55600 32.52900
C  46.35500 24.06200 32.46000
C  46.45000 23.76800 30.90100
C  47.14000 23.08900 33.32600
C  45.54800 26.46400 33.25900
C  45.96500 27.28000 34.50500
H  45.04500 27.40000 35.69500
N  46.38000 23.04600 27.85500
C  46.17200 22.14900 28.89200
C  45.97200 20.84900 28.36300
C  46.10600 20.92300 26.95000
C  46.38200 22.34000 26.69200
C  45.54700 19.69500 29.35300
C  46.04800 19.73200 25.99500
O  46.24200 19.78300 24.77100
C  45.65900 18.41400 26.50400
N  47.46200 25.07000 25.97900
C  47.20100 24.13500 25.08900
C  47.75300 24.52400 23.69900
C  47.64000 26.09100 23.86200
C  47.64600 26.26000 25.37800
C  49.21100 24.11600 23.27000
C  46.37000 26.80300 23.19000
C  45.02600 26.38600 23.62800
N  47.60700 27.03800 28.18700
C  48.00800 27.95100 27.19200
C  48.42700 29.16900 27.80900
C  48.29400 28.93400 29.16300
C  47.68500 27.64800 29.33600
C  49.11900 30.26800 27.16400
C  48.30700 29.55700 30.40200
O  48.64500 30.69400 30.74700
C  47.85400 28.48300 31.49300
C  49.08400 28.17700 32.28200
O  50.16500 27.73400 31.95500
O  48.78800 28.53400 33.53300
C  49.86000 28.29900 34.55500
H  45.83500 21.89200 31.02800
H  46.47600 22.23000 24.61700
H  48.19300 28.30200 25.21900
H  47.64500 25.71200 33.03300
H  45.30400 23.92600 32.71400
H  47.87300 22.82100 32.56600
H  46.56600 22.26800 33.75600
H  47.69200 23.45800 34.19000
H  44.70400 25.82300 33.51200
H  45.17000 27.27500 32.63600
H  46.05600 28.31300 34.16800
H  46.96700 27.02300 34.84600
H  46.22400 18.85300 29.21000
H  44.59100 19.46800 28.88100
H  45.60100 19.88500 30.42500
H  46.32300 18.09200 27.30500
H  45.94900 17.70600 25.72800
H  44.60000 18.29800 26.73600
H  46.98600 24.30500 22.95600
H  48.53000 26.46300 23.35500
H  49.97100 24.18900 24.04700
H  49.40300 24.79200 22.43600
H  49.23500 23.08400 22.91900
H  46.32700 26.67600 22.10800
H  46.28900 27.84300 23.50600
H  44.88500 26.94500 24.55200
H  44.84400 25.33800 23.86800
H  44.17200 26.64500 23.00300
H  50.15800 29.93700 27.14300
H  49.06200 31.22600 27.68000
H  48.71000 30.52400 26.18600
H  47.00000 28.80500 32.08900
H  49.50700 27.43600 35.12000
H  49.99000 29.19400 35.16300
H  50.84800 28.01700 34.19000
Mg 40.61000 41.38800 26.95600
C  39.86800 43.56800 29.58900
C  41.56600 39.09100 29.09600
C  42.01100 39.61800 24.36600
C  40.07100 44.02500 24.67100
N  40.59100 41.34900 29.06300
C  40.20000 42.29100 29.98900
C  40.16900 41.71800 31.40000
C  40.93300 40.36900 31.20600
C  40.94200 40.20500 29.72600
C  42.27000 40.16100 31.88400
C  38.74600 41.61400 32.03600
C  38.70900 41.06400 33.46600
H  37.95300 41.90400 34.48500
N  41.75400 39.53700 26.76500
C  41.90400 38.67000 27.78800
C  42.49400 37.45600 27.29600
C  42.73100 37.64100 25.92300
C  42.16600 38.98800 25.63800
C  42.83400 36.26800 28.17100
C  43.30000 36.69700 24.91500
O  43.61700 36.97700 23.72800
C  43.50600 35.23000 25.27200
N  41.00200 41.78000 24.87100
C  41.64700 40.91600 24.05300
C  41.55500 41.44600 22.52200
C  40.51800 42.61700 22.64600
C  40.44500 42.80800 24.12700
C  42.99600 41.87100 22.03700
C  39.09100 42.15200 22.17700
C  38.49200 40.88200 22.76100
N  40.00500 43.44200 27.05400
C  39.84200 44.33300 26.02900
C  39.43100 45.53300 26.59500
C  39.42300 45.36800 28.01000
C  39.78400 44.06700 28.24300
C  39.15600 46.83600 25.93100
C  39.20700 45.98800 29.33700
O  38.93200 47.12500 29.70300
C  39.39600 44.78700 30.41700
C  40.24300 45.23600 31.55400
O  41.22600 46.02900 31.45500
O  39.66700 44.79200 32.66300
C  40.07000 45.50000 33.86300
H  41.79000 38.25100 29.75700
H  42.42000 39.05300 23.52500
H  39.79000 44.77400 23.92900
H  40.67200 42.47500 32.00200
H  40.36400 39.50500 31.55100
H  42.49900 40.80900 32.73000
H  43.06200 40.23200 31.13800
H  42.35700 39.09400 32.09000
H  38.09900 41.02000 31.39000
H  38.28900 42.60300 32.04400
H  39.66400 41.20800 33.97100
H  38.39800 40.01900 33.48100
H  42.71500 36.43800 29.24100
H  43.88400 35.97700 28.18900
H  42.26200 35.37200 27.92800
H  44.31400 35.17500 26.00100
H  43.97700 34.67500 24.46100
H  42.58300 34.72100 25.55200
H  41.23200 40.63200 21.87300
H  40.99300 43.46400 22.15100
H  43.79600 41.51800 22.68800
H  43.14700 42.94900 21.97800
H  43.27300 41.44900 21.07200
H  39.29200 42.04500 21.11100
H  38.30900 42.90300 22.28500
H  39.25300 40.32600 23.30900
H  37.94200 40.24300 22.07100
H  37.70100 41.27300 23.40100
H  38.42700 46.78300 25.12200
H  40.10600 47.12900 25.48500
H  38.86200 47.61200 26.63700
H  38.38700 44.62200 30.79400
H  39.19900 45.70100 34.48800
H  40.53400 46.47500 33.71600
H  40.74300 44.82600 34.39400


