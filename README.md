compile for speed

g++ simplemlt.cpp -O3 -march=native -funroll-loops -mfpmath=sse -msse2 -ffast-math -fwhole-program -Wno-format-y2k  -fno-exceptions -fno-strict-aliasing -s -fopenmp -std=c++0x -mconsole -lole32 -luuid -lcomctl32
