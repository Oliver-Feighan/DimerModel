%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 40.54200 8.71200 29.57000
C  42.46900 10.35300 31.83700
C  38.24000 7.70200 31.93300
C  38.87500 6.96300 27.14400
C  43.10200 9.33000 27.08900
N  40.50700 8.83000 31.66000
C  41.28100 9.79100 32.39300
C  40.70700 10.01900 33.73400
C  39.62800 8.84800 33.81600
C  39.41100 8.39300 32.39000
C  40.09700 7.61800 34.74300
C  40.24300 11.48300 34.00600
C  40.64100 12.04200 35.33400
H  40.55300 11.05000 36.56000
N  38.75400 7.55700 29.53000
C  38.00600 7.28200 30.63600
C  36.89200 6.50300 30.21600
C  36.97400 6.32300 28.78000
C  38.23000 7.03300 28.36500
C  35.87900 5.90900 31.21900
C  36.08100 5.65900 27.77900
O  36.28700 5.46300 26.57600
C  34.68000 5.19900 28.33400
N  40.98100 8.10600 27.39000
C  40.09800 7.49500 26.64000
C  40.46000 7.54300 25.15700
C  41.81600 8.30500 25.16100
C  42.00400 8.59400 26.64900
C  40.56900 6.08800 24.47000
C  41.66000 9.63200 24.29500
C  40.74000 10.78700 24.76200
N  42.48900 9.50900 29.49800
C  43.34000 9.72700 28.43100
C  44.44800 10.55700 28.94600
C  44.20300 10.86700 30.25000
C  42.96700 10.18300 30.54000
C  45.69600 10.84800 28.18600
C  44.66700 11.56800 31.44700
O  45.66700 12.23600 31.65000
C  43.53700 11.23400 32.51700
C  44.17900 10.46300 33.60100
O  44.82400 9.40400 33.50900
O  43.90700 11.19500 34.73900
C  44.28800 10.53200 35.99300
H  37.46200 7.40700 32.64000
H  38.23300 6.48700 26.39900
H  43.84700 9.59300 26.33600
H  41.53500 9.82500 34.41700
H  38.70200 9.27400 34.20100
H  40.95100 7.85100 35.37900
H  40.49700 6.86600 34.06300
H  39.17400 7.30300 35.22800
H  39.15600 11.54000 33.94800
H  40.73700 12.04000 33.21000
H  39.99600 12.85500 35.66500
H  41.60800 12.52200 35.48500
H  35.14900 6.60800 31.62700
H  36.53200 5.46300 31.96900
H  35.20300 5.18600 30.76300
H  34.82100 4.32000 28.96300
H  33.98600 5.15200 27.49500
H  34.17100 5.96300 28.92200
H  39.71300 8.14400 24.63900
H  42.64100 7.73100 24.73900
H  39.67900 5.98600 23.84800
H  40.55600 5.31000 25.23400
H  41.39700 5.86900 23.79700
H  41.39900 9.47200 23.24900
H  42.67000 10.03800 24.36000
H  41.45000 11.59500 24.93900
H  40.14900 10.52600 25.63900
H  40.01100 11.10100 24.01500
H  45.97000 11.83700 28.55200
H  45.54800 10.83600 27.10600
H  46.48500 10.23100 28.61400
H  43.32700 12.24400 32.86800
H  44.09900 11.31400 36.72900
H  45.36500 10.37600 35.93100
H  43.84100 9.55300 36.16700
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


