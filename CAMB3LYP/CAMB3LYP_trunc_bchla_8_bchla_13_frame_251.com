%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5)

TDDFT excited states in gaussian

0 1
Mg 45.00000 3.44200 47.03300
C  42.74800 5.98100 46.47100
C  42.51500 1.08600 46.99500
C  47.33500 0.85800 47.41000
C  47.61400 5.55900 46.25100
N  42.89400 3.58700 46.53800
C  42.09600 4.73400 46.51200
C  40.62300 4.38600 46.28300
C  40.66000 2.87900 46.64300
C  42.13500 2.42800 46.65300
C  39.85400 2.54000 47.91600
C  40.17300 4.67400 44.83700
C  38.76200 5.33200 44.65300
H  38.03300 4.97800 43.41600
N  44.91300 1.33100 47.23000
C  43.80200 0.56100 47.28100
C  44.07100 -0.78100 47.70300
C  45.44900 -0.85000 47.96600
C  45.98000 0.51500 47.55400
C  42.99500 -1.83600 47.96800
C  46.25300 -2.09300 48.49400
O  45.64300 -3.18400 48.63300
C  47.77600 -2.03300 48.94400
N  47.15200 3.17200 46.56100
C  47.86200 2.08300 46.94400
C  49.35200 2.30000 46.87800
C  49.47200 3.78200 46.46200
C  48.00700 4.22600 46.45700
C  50.12800 1.90800 48.17900
C  50.29800 3.99400 45.10500
C  51.71400 4.54900 45.29800
N  45.19600 5.38500 46.52800
C  46.33300 6.17700 46.25500
C  45.93600 7.59000 46.19700
C  44.55800 7.53700 46.35400
C  44.18300 6.18500 46.42000
C  46.76000 8.81000 45.98900
C  43.34800 8.37000 46.40400
O  43.19200 9.59500 46.53700
C  42.11400 7.38900 46.52800
C  41.34500 7.62500 47.82100
O  41.80700 7.77900 48.95600
O  40.04700 7.69400 47.41000
C  39.07000 8.09100 48.41000
H  41.73100 0.32600 46.99100
H  48.15700 0.16300 47.59400
H  48.39400 6.28300 46.00700
H  39.99900 4.98400 46.94800
H  40.15100 2.41000 45.80100
H  39.07600 3.22300 48.25700
H  40.58100 2.53900 48.72800
H  39.36600 1.56800 47.84800
H  40.24600 3.80100 44.18900
H  40.97000 5.28600 44.41400
H  38.87100 6.41600 44.68600
H  38.11500 4.99000 45.46200
H  43.14300 -2.15000 49.00100
H  42.99500 -2.57100 47.16300
H  42.06000 -1.27600 47.95400
H  48.36300 -2.15900 48.03500
H  47.92900 -2.88500 49.60600
H  47.87200 -1.14900 49.57500
H  49.74800 1.66300 46.08700
H  49.86100 4.46700 47.21500
H  50.84800 2.67500 48.46300
H  50.67000 0.96800 48.07700
H  49.41000 1.70300 48.97200
H  49.82900 4.72000 44.44100
H  50.30900 3.07000 44.52700
H  52.39600 4.06600 44.59800
H  52.05000 4.35300 46.31600
H  51.79100 5.62600 45.15000
H  47.80500 8.50000 45.98500
H  46.42200 9.60500 46.65400
H  46.58100 9.11700 44.95800
H  41.44200 7.45600 45.67200
H  39.19800 7.59700 49.37400
H  38.01800 7.98300 48.14600
H  39.14800 9.16300 48.59300
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


