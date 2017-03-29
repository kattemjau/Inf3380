CC=gcc
CFLAGS= -Wall -Wextra -std=c99 -g -D_XOPEN_SOURCE

NAME=Oblig1

All:compile run

compile: Oblig1.c
		$(CC) $(CFLAGS) $< -o $(NAME)

run:
	./$(NAME) $(HOST) $(PORT)

clean:
	rm -f $(NAME)

check: compile
	valgrind --leak-check=full ./$(NAME) $(HOST) $(PORT)
