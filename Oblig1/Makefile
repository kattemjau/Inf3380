CC = mpicc
CFLAGS = -Wall -Wextra -std=c99 -O2
LDFLAGS = -L ./simple-jpeg/ -l simplejpeg

PROJ = serial_main
# parallel_main
PARR = parallel_main

KAPPA=0.2
ITERS=100
INPUTF=mona_lisa_noisy.jpg
OUTPUTF=out.jpg


all: clean compile run_s

compile: clean simple-jpeg serial_main parallel

serial_main: $(PROJ)
	$(CC) $(CFLAGS) serial_main.c -o $@ $(LDFLAGS)

run_s:
	./$(PROJ) $(KAPPA) $(ITERS) $(INPUTF) $(OUTPUTF)
run_p:
	./$(PARR) $(KAPPA) $(ITERS) $(INPUTF) $(OUTPUTF)

parallel: $(PARR)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

simple-jpeg :
	$(MAKE) -C /simple-jpeg

obj-clean:
	$(RM) *.o

exec-clean:
	$(RM) $(PROJ)

autosave-clean:
	$(RM) *~

clean:
	$(MAKE) obj-clean
	$(MAKE) exec-clean
	$(MAKE) autosave-clean
