%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

TDDFT excited states in gaussian

0 1
Mg 46.26900 25.13700 28.87000
C  47.26400 27.26100 31.45200
C  45.06800 22.98000 31.10200
C  46.59300 22.69600 26.49900
C  47.58600 27.40300 26.57500
N  46.02800 25.20400 31.05400
C  46.48600 26.15500 31.85900
C  45.98900 25.86800 33.32200
C  45.39300 24.44600 33.19000
C  45.51200 24.15500 31.69700
C  46.24000 23.30800 33.94900
C  44.92100 26.90700 33.89400
C  45.18000 27.25800 35.31900
H  44.05100 27.15600 36.33300
N  45.85600 23.12600 28.82000
C  45.33200 22.40200 29.85800
C  44.87200 21.10600 29.34000
C  45.32500 21.03600 27.99300
C  46.00600 22.31800 27.74700
C  44.15900 20.03600 30.22700
C  45.21800 19.84800 26.99000
O  45.63600 19.88000 25.83800
C  44.54500 18.61000 27.44900
N  46.98900 25.03900 26.90500
C  47.10000 23.90200 26.09900
C  47.90400 24.25000 24.81900
C  47.61600 25.80400 24.63200
C  47.36900 26.07600 26.14300
C  49.38900 23.83600 24.90100
C  46.42100 26.18700 23.80600
C  45.03400 26.12100 24.37300
N  47.05200 27.04800 28.89500
C  47.58500 27.84200 27.90600
C  48.07800 29.11700 28.43200
C  48.07000 28.90600 29.88400
C  47.43900 27.61600 30.07900
C  48.67900 30.18600 27.67300
C  48.43700 29.46600 31.17900
O  49.02500 30.43400 31.48700
C  47.84100 28.42700 32.26500
C  48.81200 28.03800 33.35100
O  49.55700 27.06300 33.31900
O  48.64800 28.91100 34.42600
C  49.48500 28.68900 35.61900
H  44.47200 22.24900 31.65300
H  46.68200 21.86300 25.79800
H  47.89200 28.10500 25.79700
H  46.81500 25.79600 34.03000
H  44.37900 24.31000 33.56600
H  45.81000 23.21000 34.94600
H  47.29000 23.54400 34.12000
H  46.10500 22.37300 33.40400
H  43.90600 26.51300 33.84000
H  44.83700 27.84500 33.34500
H  45.39600 28.32500 35.27200
H  46.03100 26.79100 35.81400
H  44.95600 19.41500 30.63500
H  43.40600 19.48400 29.66400
H  43.77200 20.42600 31.16900
H  44.70200 18.21100 28.45100
H  44.79100 17.78900 26.77600
H  43.48800 18.66000 27.18600
H  47.38400 23.69600 24.03800
H  48.47400 26.33800 24.22400
H  50.12500 24.61500 24.70200
H  49.61300 22.90600 24.37800
H  49.63200 23.53200 25.91900
H  46.32100 25.62700 22.87700
H  46.56100 27.21900 23.48400
H  44.87400 27.14600 24.70700
H  44.90600 25.48200 25.24700
H  44.33600 25.84200 23.58300
H  48.36400 31.10800 28.16200
H  48.25000 30.14800 26.67200
H  49.76500 30.10200 27.69900
H  47.00200 28.92500 32.75100
H  50.51900 28.64400 35.27600
H  49.25400 27.73000 36.08200
H  49.44300 29.39800 36.44600
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


