ARGS= -Wall -Wextra -std=c99 -O2 -g
FILES=small_matrix_a.bin small_matrix_b.bin c.bin
# FILES2=mina.bin minb.bin minc.bin

run: clear compile
	mpirun -np 5 ./oblig $(FILES)

compile:
	mpicc Oblig2.c -o oblig $(ARGS)

clear:
	clear

clean:
	rm -r oblig
	rm -r c.bin

check: compile
	valgrind --leak-check=full ./oblig $(FILES)

	#checking purposes only
diff: run
	vimdiff small_matrix_c.bin c.bin

git: clean
	git add .
	git commit -m "something is working"
	git push
