%nproc=24
%mem=175gb
#p PBE1PBE/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 40.74100 8.40300 29.69900
C  42.66000 10.02500 32.01600
C  38.73200 7.21500 32.07700
C  39.17900 6.75200 27.36300
C  43.28300 9.22800 27.34100
N  40.79300 8.45600 31.86800
C  41.50500 9.39200 32.54400
C  40.94800 9.63500 33.92000
C  39.85300 8.49700 34.00200
C  39.82700 7.94800 32.59200
C  40.10000 7.44100 35.10600
C  40.44100 11.18100 34.18500
C  40.87000 11.88200 35.48300
H  41.23900 10.96600 36.69200
N  39.13800 7.34400 29.69200
C  38.39900 6.91700 30.74900
C  37.34700 6.00300 30.37600
C  37.33700 6.02300 29.01400
C  38.65500 6.62800 28.64700
C  36.61000 5.29800 31.50000
C  36.33400 5.36600 28.04500
O  36.44200 5.38700 26.84100
C  35.20700 4.65100 28.73000
N  41.30000 7.91900 27.67200
C  40.36700 7.31400 26.91500
C  40.80200 7.26400 25.44700
C  41.95100 8.32700 25.38700
C  42.23200 8.49800 26.89400
C  41.21600 5.86100 24.82900
C  41.58000 9.72400 24.76100
C  40.48400 10.59900 25.31000
N  42.57200 9.48100 29.57000
C  43.55300 9.66900 28.59100
C  44.67900 10.46300 29.09900
C  44.36400 10.63700 30.45500
C  43.08700 10.02200 30.64800
C  45.85600 11.01300 28.33500
C  44.78800 11.17000 31.66800
O  45.80500 11.75000 32.01500
C  43.61600 11.01200 32.67100
C  44.29700 10.49000 33.87800
O  44.63700 9.31800 34.05700
O  44.61300 11.49900 34.74100
C  45.54500 11.24000 35.84300
H  37.99000 7.00700 32.85100
H  38.68100 6.27000 26.52000
H  43.95700 9.44200 26.50800
H  41.68000 9.37000 34.68300
H  38.87700 8.95900 34.14900
H  40.48500 6.53200 34.64400
H  39.14600 7.20600 35.57800
H  40.83000 7.87400 35.79000
H  39.35300 11.12200 34.15300
H  40.80700 11.84200 33.40000
H  40.04300 12.46600 35.88600
H  41.65400 12.60800 35.26800
H  37.22900 4.66200 32.13300
H  35.89000 4.51000 31.28000
H  36.12800 6.14400 31.99000
H  34.67500 5.32500 29.40100
H  35.72300 3.80100 29.17600
H  34.37800 4.27900 28.12900
H  39.88500 7.56300 24.93800
H  42.83200 7.86300 24.94300
H  41.11200 5.08600 25.58800
H  42.28100 5.83900 24.60100
H  40.70000 5.56000 23.91600
H  41.36300 9.40900 23.74000
H  42.43100 10.40400 24.73900
H  39.92000 10.04700 26.06300
H  39.77700 10.86400 24.52400
H  40.79900 11.50000 25.83600
H  46.78200 10.48700 28.56300
H  45.97200 12.07500 28.55200
H  45.60900 10.77100 27.30100
H  43.07400 11.94700 32.81400
H  46.05900 12.14500 36.16900
H  46.18500 10.38200 35.64300
H  45.04400 10.86200 36.73400
Mg 26.06500 50.95300 26.70600
C  23.95100 52.21500 29.08800
C  28.12100 49.83900 29.10000
C  28.15700 49.98600 24.24600
C  24.06600 52.63600 24.17800
N  26.14100 51.27800 28.89700
C  25.13500 51.71000 29.66700
C  25.49500 51.52800 31.23000
C  26.85400 50.74300 31.10400
C  27.11600 50.62500 29.58900
C  28.00100 51.35300 31.93300
C  24.40800 50.76600 32.04000
C  24.78700 50.21400 33.45300
H  24.26800 51.00700 34.71200
N  27.93000 50.07200 26.68200
C  28.56100 49.60000 27.73700
C  29.78100 48.94300 27.29300
C  29.90300 49.05100 25.87500
C  28.62800 49.72200 25.53900
C  30.86100 48.38000 28.20200
C  31.03700 48.64600 24.93900
O  31.09500 49.07700 23.79100
C  32.08400 47.71600 25.44400
N  26.12000 51.31300 24.51700
C  27.05500 50.76500 23.75000
C  26.75200 50.94300 22.25300
C  25.36900 51.58600 22.31600
C  25.20000 51.88900 23.75800
C  27.77300 51.70600 21.39600
C  24.29900 50.82000 21.58000
C  23.63800 51.55200 20.35800
N  24.37700 52.10200 26.60600
C  23.61600 52.65500 25.52200
C  22.39000 53.25700 25.99800
C  22.52400 53.06300 27.40100
C  23.71600 52.38900 27.69500
C  21.20200 53.78500 25.27400
C  21.72000 53.18300 28.59100
O  20.53300 53.46900 28.80800
C  22.59000 52.78500 29.75700
C  22.95900 54.06000 30.54700
O  23.29600 55.16300 30.09000
O  22.85900 53.76400 31.87400
C  23.13500 54.86100 32.85500
H  28.80700 49.32300 29.77500
H  28.85800 49.89200 23.41400
H  23.44600 52.99600 23.35500
H  25.65400 52.54800 31.57800
H  26.67500 49.72400 31.44700
H  28.52200 50.52300 32.40900
H  27.63000 51.99000 32.73600
H  28.68100 51.92100 31.29800
H  24.23600 49.84500 31.48300
H  23.53400 51.40200 32.18000
H  25.86900 50.31100 33.54800
H  24.40200 49.19900 33.54500
H  30.59600 48.49500 29.25300
H  31.79900 48.90500 28.02300
H  31.04200 47.32300 28.01100
H  31.53500 46.89000 25.89700
H  32.83600 48.00900 26.17600
H  32.72300 47.32900 24.65100
H  26.66800 49.94900 21.81500
H  25.43600 52.54200 21.79700
H  27.79000 51.15500 20.45600
H  28.73100 51.77600 21.91100
H  27.54500 52.76300 21.25900
H  23.50000 50.72800 22.31500
H  24.63300 49.85300 21.20500
H  24.06100 52.55100 20.24900
H  22.57600 51.65700 20.58200
H  23.65600 51.00200 19.41700
H  21.27300 53.54700 24.21300
H  21.29800 54.86900 25.32400
H  20.34800 53.31100 25.76000
H  21.96100 52.11200 30.34000
H  24.03500 55.44400 32.65800
H  23.15200 54.54000 33.89600
H  22.32000 55.57300 32.72600


