CC = mpicc
GCC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2 -g
LDFLAGS = -L ./simple-jpeg/ -l simplejpeg

PROJ = serial_main
# parallel_main
PARR = parallel_main

KAPPA=0.2
ITERS=100

INPUTF=mona_lisa_noisy.jpg
OUTPUTF=out.jpg
PROCS=8

all: parallel run_p

compile: clean simple-jpeg serial_main parallel

serial_main:
	$(GCC) $(CFLAGS) serial_main.c -o $@ $(LDFLAGS)

run_s:
	./$(PROJ) $(KAPPA) $(ITERS) $(INPUTF) $(OUTPUTF)
run_p:
	mpirun -np $(PROCS) ./parallel_main $(KAPPA) $(ITERS) $(INPUTF) $(OUTPUTF)

parallel:
	$(CC) $(CFLAGS) parallel_main.c -o parallel_main $(LDFLAGS)

simple-jpeg :
	$(MAKE) -C /simple-jpeg

obj-clean:
	$(RM) *.o

exec-clean:
	$(RM) serial_main
	$(RM) parallel_main

autosave-clean:
	$(RM) *~

clean:
	$(MAKE) obj-clean
	$(MAKE) exec-clean
	$(MAKE) autosave-clean

check: parallel
		mpirun -np $(PROCS) valgrind --leak-check=full ./parallel_main $(KAPPA) $(ITERS) $(INPUTF) $(OUTPUTF)

check_s: serial_main
		valgrind --leak-check=full ./serial_main $(KAPPA) $(ITERS) $(INPUTF) $(OUTPUTF)
