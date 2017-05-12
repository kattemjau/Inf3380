ARGS= -Wall -Wextra -std=c99 -O2 -g
FILES=small_matrix_a.bin small_matrix_b.bin small_matrix_c.bin

run: clear compile
	./hei $(FILES)

compile:
	gcc Oblig2.c -o hei $(ARGS)

clear:
	clear

clean:
	rm -r hei

check:
	valgrind --leak-check=full ./hei $(FILES)
