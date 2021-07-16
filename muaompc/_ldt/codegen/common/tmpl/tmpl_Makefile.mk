CC = gcc
FLAGS = -g -Wall -Wstrict-prototypes -pedantic
OPT = -O0 -funroll-loops
STD = -std=c89

OBJ := $(patsubst %.c,%.o,$(wildcard *.c))
INC := ./include

main: $(OBJ)
	$(CC) $(FLAGS) $(OPT) $(STD) -I$(INC) $(OBJ) -lm -o main

lib{prefix}: $(OBJ)
	ar rcs lib{prefix}.a $(OBJ)

$(OBJ): %.o: %.c
	$(CC) $(FLAGS) $(OPT) $(STD) -I$(INC) -c $< -o $@

.PHONY: clean

clean:
	rm main $(OBJ) lib{prefix}.a
