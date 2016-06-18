compile for speed

g++ simplemlt.cpp -O3 -march=native -funroll-loops -mfpmath=sse -msse2 -msse3 -mssse3  -ffast-math -fwhole-program -Wno-format-y2k  -fno-exceptions -fno-strict-aliasing -s -fopenmp -std=c++0x -mconsole -lole32 -luuid -lcomctl32

try run with command
mlt.exe -mutations 3000000 -width 600 -height 480

compared with original version, new version do test command in 54s, original version - 65s. result are same
