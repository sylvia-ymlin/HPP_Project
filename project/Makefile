# Build executable
galsim: galsim.c barnes_hut.c kmeans.c
	gcc -Wall -Wextra -O3 -march=native -ffast-math -fopenmp galsim.c barnes_hut.c kmeans.c -lm -o galsim

# Clean target
clean:
	rm -f galsim result.gal
