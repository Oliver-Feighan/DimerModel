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
Mg 34.81300 49.43200 25.14500
C  35.01600 47.65700 28.22400
C  33.32400 52.02800 26.83800
C  34.65200 50.90000 22.18900
C  35.90600 46.46600 23.43400
N  34.14900 49.79800 27.34500
C  34.46100 48.93800 28.38200
C  33.98800 49.65100 29.69500
C  33.38500 51.08500 29.24000
C  33.61700 50.95700 27.73100
C  31.94200 51.39000 29.69300
C  35.01900 49.70600 30.89700
C  34.55900 49.07300 32.17800
H  33.11300 48.96700 32.58500
N  34.17200 51.35300 24.58700
C  33.58900 52.24400 25.43200
C  33.45200 53.46600 24.68200
C  33.72500 53.14400 23.33500
C  34.19500 51.72200 23.26900
C  32.90900 54.73200 25.36600
C  33.46900 53.99900 22.13100
O  33.64100 53.56300 20.96800
C  32.94500 55.37500 22.31400
N  35.18300 48.68800 23.08400
C  35.11600 49.51200 22.09700
C  35.57200 48.95700 20.64400
C  35.95000 47.52900 21.06600
C  35.68700 47.53300 22.55200
C  34.36900 49.06600 19.63700
C  37.40800 47.05400 20.73300
C  38.54700 47.67100 21.57900
N  35.21400 47.49800 25.67400
C  35.61800 46.42200 24.84900
C  35.92400 45.27900 25.72000
C  35.71000 45.72300 27.04700
C  35.31700 47.10200 26.97200
C  36.42700 43.88500 25.23900
C  35.78900 45.30100 28.34400
O  36.05400 44.22400 28.83800
C  35.40700 46.52600 29.22100
C  34.34200 46.13400 30.11500
O  33.13700 46.19100 29.84500
O  34.87800 45.46700 31.22300
C  33.86400 44.85200 32.11500
H  32.63300 52.79500 27.19500
H  34.43900 51.33900 21.21200
H  36.24500 45.54000 22.96500
H  33.26200 48.90600 30.01800
H  33.96500 51.87500 29.71500
H  31.90300 52.32700 30.25000
H  31.43700 50.62000 30.27600
H  31.39300 51.49300 28.75700
H  35.30500 50.73900 31.09600
H  35.87800 49.12400 30.56300
H  34.83800 49.76200 32.97600
H  35.03100 48.09400 32.25300
H  32.96100 54.66600 26.45300
H  31.87600 54.89800 25.06300
H  33.55400 55.53200 25.00400
H  33.66100 55.90400 22.94300
H  31.99400 55.37100 22.84800
H  32.75700 55.84600 21.35000
H  36.44300 49.49500 20.26900
H  35.19100 46.83200 20.71000
H  33.44700 48.97400 20.21100
H  34.41800 48.21400 18.95800
H  34.29100 49.96700 19.02800
H  37.57700 47.10800 19.65800
H  37.43700 45.97400 20.87400
H  38.40300 48.69400 21.92500
H  39.48100 47.54600 21.03000
H  38.64300 47.01100 22.44100
H  37.30600 43.53600 25.78100
H  36.62200 43.84300 24.16800
H  35.74400 43.11800 25.60500
H  36.33300 46.69700 29.76900
H  33.38300 44.00700 31.62300
H  33.04700 45.51600 32.39500
H  34.36200 44.67800 33.06900


