%nproc=24
%mem=175gb
#p cam-b3lyp/Def2SVP td=(nstates=5) density=(transition=1)

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
Mg -5.36100 24.76000 27.02500
C  -3.64200 26.63800 29.43600
C  -6.40000 22.65900 29.54700
C  -6.94500 22.99100 24.73300
C  -4.07500 26.91800 24.56900
N  -5.09300 24.76200 29.22400
C  -4.33400 25.60400 30.03700
C  -4.64200 25.43000 31.57600
C  -5.37000 24.10100 31.51700
C  -5.66000 23.80300 29.98800
C  -4.43700 22.92100 32.05400
C  -5.65800 26.61300 31.90200
C  -5.45100 27.38800 33.27700
H  -4.92400 26.43700 34.32100
N  -6.35800 23.02700 27.11800
C  -6.68000 22.31300 28.19300
C  -7.48400 21.16300 27.87400
C  -7.81800 21.31800 26.51400
C  -6.97100 22.47200 26.02900
C  -7.92400 20.11100 28.82500
C  -8.85200 20.60400 25.66800
O  -9.21200 21.01500 24.55900
C  -9.41800 19.25800 26.12900
N  -5.43800 24.92500 24.91100
C  -6.26300 24.05500 24.23800
C  -6.35200 24.50400 22.78100
C  -5.50900 25.87000 22.72900
C  -5.01000 25.97700 24.14600
C  -5.90500 23.33400 21.83500
C  -6.23300 27.10700 22.22100
C  -5.47900 27.85000 21.07400
N  -4.10000 26.41500 26.99800
C  -3.63700 27.19900 25.89500
C  -2.75500 28.21500 26.34700
C  -2.74300 28.01700 27.73800
C  -3.53100 26.89000 28.05100
C  -2.03600 29.28900 25.59600
C  -2.22900 28.50400 28.97800
O  -1.46700 29.40000 29.25600
C  -2.89200 27.68300 30.12000
C  -1.83900 27.24500 31.07200
O  -0.96900 26.39800 30.89000
O  -2.04600 27.91300 32.21100
C  -1.18400 27.58800 33.32800
H  -6.90800 22.00000 30.25400
H  -7.48200 22.35700 24.02300
H  -3.70300 27.64900 23.84800
H  -3.68500 25.58900 32.07200
H  -6.31300 24.16400 32.05900
H  -3.44500 23.27200 32.33900
H  -4.33700 22.17800 31.26300
H  -4.99100 22.51200 32.89900
H  -6.69800 26.29000 31.93800
H  -5.59700 27.30300 31.06000
H  -6.41300 27.65600 33.71400
H  -4.98200 28.36000 33.12200
H  -9.00400 20.22200 28.93100
H  -7.47900 20.25200 29.80900
H  -7.70000 19.10000 28.48600
H  -9.79600 18.88600 25.17700
H  -10.18200 19.34000 26.90300
H  -8.70400 18.49400 26.43600
H  -7.35300 24.74000 22.42000
H  -4.63900 25.71000 22.09200
H  -6.65200 23.29700 21.04200
H  -5.78500 22.37600 22.34100
H  -5.04100 23.62000 21.23400
H  -6.46800 27.83800 22.99400
H  -7.10900 26.70500 21.71000
H  -4.60800 27.27900 20.75200
H  -5.00600 28.76500 21.43000
H  -6.13000 28.01900 20.21600
H  -2.85300 29.97800 25.38000
H  -1.55300 28.80400 24.74800
H  -1.32200 29.82700 26.21800
H  -3.47700 28.43300 30.65400
H  -0.13900 27.68100 33.03100
H  -1.23800 26.56800 33.71000
H  -1.41700 28.31200 34.10900


