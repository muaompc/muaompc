CC = gcc
FLAGS = -g -Wall -Wstrict-prototypes -pedantic
OPT = -O0 -funroll-loops
STD = -std=c89

OBJ := $(patsubst %.c,%.o,$(wildcard *.c))
INC := ./include

libstocmpc: $(OBJ)
	ar rcs libstocmpc.a $(OBJ)

$(OBJ): %.o: %.c
	$(CC) $(FLAGS) $(OPT) $(STD) -I$(INC) -c $< -o $@

.PHONY: clean

clean:
	rm $(OBJ) libstocmpc.a
