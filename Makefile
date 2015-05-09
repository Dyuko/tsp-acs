CC = gcc

SRC   = src/tsp-acs.c
FLAGS = -lm -lpthread -Wall
OUT   = tsp-acs

compile:
	$(CC) $(FLAGS) $(SRC) -o $(OUT)
	
test:
	time ./$(OUT) data/a280.tsp 10 0 
	time ./$(OUT) data/a280.tsp 10 1