CC = gcc
FLAGS = -g -Wall -Wstrict-prototypes -pedantic
OPT = -O0 -funroll-loops
STD = -std=c89

OBJ := $(patsubst %.c,%.o,$(wildcard *.c)) $(patsubst %.c,%.o,$(wildcard ../../src/*.c))
INC := ../../src/include
INCD := .

main: $(OBJ)
	$(CC) $(FLAGS) $(OPT) $(STD) -I$(INC) -I$(INCD) $(OBJ) -lm -o main

$(OBJ): %.o: %.c
	$(CC) $(FLAGS) $(OPT) $(STD) -I$(INC) -c $< -o $@

.PHONY: clean

clean:
	rm main
	rm $(OBJ)
