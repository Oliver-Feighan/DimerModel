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
Mg 7.14200 57.54600 41.51700
C  6.28900 54.20700 41.28700
C  10.28800 56.72100 40.38900
C  7.81200 60.89700 41.16900
C  3.94100 58.30400 42.53200
N  8.10000 55.73200 40.72500
C  7.55000 54.46800 40.80800
C  8.51900 53.39300 40.32300
C  9.88500 54.17400 40.37000
C  9.39700 55.61300 40.45200
C  10.94100 53.73800 41.48400
C  8.18100 52.71900 38.95500
C  8.31100 51.21100 38.91500
H  9.48100 50.75900 38.05500
N  8.80400 58.60500 40.85800
C  9.98200 58.10200 40.47200
C  10.93200 59.16600 40.19700
C  10.33500 60.40200 40.38100
C  8.90800 60.02900 40.77200
C  12.41400 58.88100 39.65800
C  11.01900 61.80700 40.27800
O  12.22600 61.86700 39.90900
C  10.28000 63.06100 40.55000
N  5.92900 59.42300 41.58100
C  6.46000 60.62400 41.57700
C  5.53200 61.67900 42.14200
C  4.20400 60.94400 42.52900
C  4.70000 59.43900 42.28400
C  6.08200 62.75400 43.13300
C  2.88100 61.24800 41.76800
C  1.64800 60.78800 42.51800
N  5.41300 56.46200 41.84200
C  4.16700 56.89800 42.34900
C  3.34900 55.73600 42.56600
C  4.13900 54.62000 42.19700
C  5.33900 55.15600 41.73800
C  2.01600 55.69000 43.29800
C  4.11700 53.24300 42.14700
O  3.28800 52.39200 42.32200
C  5.50100 52.88200 41.39400
C  5.07100 52.59100 39.95700
O  4.83200 53.42500 39.10500
O  4.99200 51.21300 39.83600
C  4.48800 50.71300 38.56100
H  11.34700 56.50800 40.23000
H  8.09700 61.94400 41.04900
H  2.94700 58.37800 42.97800
H  8.57300 52.63000 41.09900
H  10.29300 54.14500 39.35900
H  11.53300 54.64700 41.59000
H  11.67500 53.00700 41.14400
H  10.53800 53.39000 42.43500
H  8.89800 53.17500 38.27200
H  7.26400 53.03800 38.46000
H  7.37100 50.81600 38.52900
H  8.48000 50.78000 39.90200
H  12.65000 57.83100 39.48000
H  13.15600 59.27900 40.35100
H  12.59000 59.28000 38.66000
H  9.90800 63.04900 41.57400
H  9.55300 63.27700 39.76700
H  11.06400 63.81300 40.63000
H  5.09100 62.35700 41.41100
H  4.06100 60.88100 43.60800
H  7.16900 62.67900 43.08000
H  5.65900 62.68700 44.13600
H  5.85900 63.78200 42.84700
H  2.99200 60.79700 40.78200
H  2.91000 62.32500 41.60300
H  0.83700 61.47100 42.26400
H  1.73700 60.78400 43.60500
H  1.32700 59.82700 42.11600
H  1.67100 54.68100 43.52400
H  1.28900 56.22800 42.69100
H  2.06800 56.22200 44.24800
H  6.04900 52.10100 41.92200
H  3.44600 51.01000 38.44100
H  4.61600 49.63500 38.46900
H  5.02300 51.19800 37.74500


