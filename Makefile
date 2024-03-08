
# Build executable
galsim: galsim.c
	gcc -Wall -Wextra -O3 -march=native -ffast-math -fopenmp galsim.c -lm -o galsim

# Clean target
clean:
	rm -f galsim result.gal
