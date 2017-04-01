CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2
LDFLAGS = -L ./simple-jpeg/ -l simplejpeg

PROJ = serial_main
# parallel_main
PARR = parallel_main.o
SERR = serial_main.o

KAPPA=0.2
ITERS=1
INPUTF=mona_lisa_noisy.jpg
OUTPUTF=out.jpg


all: simple-jpeg $(PROJ) run_s

serial_main : $(SERR)
	gcc $(CFLAGS) $< -o $@ $(LDFLAGS)

run_s:
	./$(SERR) $(KAPPA) $(ITERS) $(INPUTF) $(OUTPUTF)
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
