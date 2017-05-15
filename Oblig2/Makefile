ARGS= -Wall -Wextra -std=c99 -O2 -g
FILES=small_matrix_a.bin small_matrix_b.bin c.bin
# FILES2=mina.bin minb.bin minc.bin

run: clear compile
	./hei $(FILES)

compile:
	gcc Oblig2.c -o hei $(ARGS)

clear:
	clear

clean:
	rm -r hei
	rm -r c.bin

check:
	valgrind --leak-check=full ./hei $(FILES)
