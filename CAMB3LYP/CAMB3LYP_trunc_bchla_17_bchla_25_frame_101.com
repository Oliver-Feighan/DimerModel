%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 29.46400 59.14500 41.03800
C  26.57800 57.38900 40.28900
C  31.34900 56.53100 39.78200
C  32.15000 61.18900 41.08600
C  27.35400 61.82000 41.94900
N  28.99200 57.29200 39.87100
C  27.77200 56.74800 39.89300
C  27.77900 55.38900 39.12700
C  29.32900 54.99700 39.21800
C  29.96700 56.34400 39.73200
C  29.58100 53.72700 40.06600
C  27.37400 55.59400 37.68500
C  26.46000 54.52700 37.10300
H  26.57800 54.02100 35.61600
N  31.53300 58.78800 40.65000
C  32.11000 57.71400 40.08400
C  33.51500 57.95800 40.03300
C  33.80600 59.24700 40.57400
C  32.47900 59.81600 40.84000
C  34.47700 56.89100 39.66100
C  35.16200 59.85800 40.59300
O  36.13000 59.21400 40.22400
C  35.39700 61.26400 41.13200
N  29.71300 61.28200 41.56500
C  30.94400 61.84600 41.41400
C  30.85900 63.33400 41.52000
C  29.29200 63.52200 41.79300
C  28.69900 62.13600 41.84900
C  31.79900 63.96200 42.56600
C  28.51800 64.48000 40.78300
C  27.64200 65.63900 41.31400
N  27.43700 59.54700 41.14000
C  26.74400 60.62700 41.58900
C  25.41200 60.28000 41.81900
C  25.25400 58.98400 41.34800
C  26.48900 58.61100 40.88000
C  24.32700 61.25500 42.36300
C  24.37200 57.92800 41.16800
O  23.15900 57.91800 41.34600
C  25.14900 56.75500 40.39700
C  25.25800 55.58200 41.18700
O  25.99300 55.31600 42.14100
O  24.33400 54.68200 40.66500
C  24.17500 53.35900 41.24400
H  31.95200 55.65500 39.53500
H  32.91900 61.96200 41.15300
H  26.80500 62.68200 42.33400
H  27.08300 54.71900 39.63100
H  29.81700 54.88800 38.25000
H  28.66900 53.37400 40.54700
H  30.34000 53.93800 40.82000
H  30.00600 52.97200 39.40500
H  28.16800 55.74800 36.95400
H  26.92700 56.58800 37.67000
H  25.44000 54.85100 37.30900
H  26.66400 53.69100 37.77200
H  34.81300 56.52500 40.63100
H  35.29700 57.28300 39.05900
H  34.02200 56.08300 39.08800
H  35.10600 61.93900 40.32700
H  36.47900 61.37300 41.20200
H  34.84400 61.42100 42.05800
H  31.14300 63.61000 40.50500
H  29.10900 64.04800 42.73000
H  32.59600 64.44300 41.99800
H  32.19100 63.23400 43.27700
H  31.26900 64.75100 43.09900
H  27.85100 63.86700 40.17700
H  29.18800 65.01800 40.11200
H  26.69000 65.11600 41.40400
H  27.54400 66.41700 40.55800
H  27.77900 66.00900 42.33000
H  24.06200 61.95200 41.56900
H  24.63400 61.90700 43.18100
H  23.42900 60.71200 42.65900
H  24.61600 56.61800 39.45600
H  23.13300 53.07400 41.10200
H  24.52900 53.25500 42.27000
H  24.72400 52.65900 40.61400
Mg -2.59000 34.01700 26.79400
C  -3.62000 32.41400 29.77300
C  -1.15300 36.45800 28.71900
C  -2.20600 35.82500 24.04600
C  -4.51700 31.70900 25.00200
N  -2.55200 34.47200 29.00300
C  -2.95200 33.59300 30.01300
C  -2.37700 34.05700 31.37400
C  -1.60200 35.37500 31.00700
C  -1.71600 35.46600 29.46000
C  -1.91400 36.61400 31.82800
C  -1.50900 33.04200 32.13400
C  -1.98900 32.76400 33.58500
H  -0.95800 32.99700 34.66700
N  -1.70200 35.86000 26.44900
C  -1.19200 36.74900 27.34900
C  -0.55100 37.80700 26.60500
C  -0.98500 37.71000 25.21400
C  -1.68000 36.43200 25.23200
C  0.22600 39.02700 27.32500
C  -1.03600 38.70200 24.13300
O  -1.40600 38.53400 23.00600
C  -0.57600 40.11000 24.49900
N  -3.46600 33.87800 24.82300
C  -3.08300 34.73100 23.85300
C  -3.52500 34.22200 22.47900
C  -3.92700 32.76700 22.76900
C  -3.94500 32.73700 24.28100
C  -4.74000 35.08300 21.90100
C  -2.77000 31.79500 22.20500
C  -1.44300 31.78800 23.01200
N  -3.75200 32.29300 27.19100
C  -4.44100 31.48200 26.37200
C  -5.05500 30.45500 27.17000
C  -4.72100 30.75200 28.47800
C  -3.96900 31.93000 28.47700
C  -5.92100 29.29000 26.73000
C  -4.97000 30.27500 29.84100
O  -5.63000 29.44100 30.39400
C  -4.20300 31.28400 30.73600
C  -5.22200 31.86200 31.64700
O  -6.22400 32.38000 31.22200
O  -4.94100 31.66200 32.98900
C  -6.05600 31.94700 33.95300
H  -0.52600 37.10900 29.33200
H  -2.00600 36.32000 23.09300
H  -4.95700 30.83400 24.51900
H  -3.20300 34.37600 32.01000
H  -0.55500 35.23200 31.27500
H  -2.84600 36.49700 32.37900
H  -1.95200 37.53200 31.24100
H  -1.13500 36.66900 32.58800
H  -0.45800 33.28400 32.29200
H  -1.47400 32.09100 31.60300
H  -2.46500 31.78700 33.66800
H  -2.77300 33.48500 33.81300
H  1.18600 39.15100 26.82400
H  0.47900 38.88500 28.37500
H  -0.26700 39.99700 27.38900
H  0.46900 39.96700 24.77500
H  -1.19700 40.39100 25.34900
H  -0.59300 40.82400 23.67500
H  -2.71600 34.17900 21.74900
H  -4.90200 32.39200 22.45700
H  -4.66500 35.04300 20.81400
H  -4.79700 36.14200 22.15500
H  -5.62600 34.53900 22.22800
H  -2.53100 31.98600 21.15900
H  -3.07000 30.75600 22.34200
H  -0.61000 31.26900 22.53900
H  -1.63900 31.28300 23.95800
H  -1.17900 32.82600 23.21700
H  -5.41600 28.57400 26.08200
H  -6.72000 29.76700 26.16300
H  -6.40800 28.78800 27.56600
H  -3.52700 30.64600 31.30500
H  -5.48500 32.48600 34.70900
H  -6.48400 31.04500 34.39000
H  -6.92400 32.52400 33.63500


