MPC = ./myprb_mydat_mpc/
MPC_SRC := $(wildcard $(MPC)*.c)
MPC_HDR := $(wildcard $(MPC)*.h)
MPC_OBJ := $(wildcard $(MPC)*.o)

CC = gcc
ASAN = -fsanitize=address  # print more information when stack smashing occurs

CFLAGS = -g $(ASAN)
LDFLAGS = $(ASAN) 


mpcmydatmain_staticmem: mpcmydatmain_staticmem.o $(MPC)libmuaompc.a 
	$(CC) mpcmydatmain_staticmem.o $(LDFLAGS) -L$(MPC) -lmuaompc -lm -o $@

mpcmydatmain_staticmem.o: mpcmydatmain_staticmem.c
	$(CC) $(CFLAGS) -c $< -o $@ -I$(MPC)

$(MPC)%.o: $(MPC)%.c $(MPC)%.h
	$(CC) $(CFLAGS) -c $< -o $@
	@echo "target: $@ prereq: $<"

$(MPC)libmuaompc.a: $(MPC_SRC) $(MPC_HDR)
	cd $(MPC); $(CC) $(CFLAGS) -c *.c; ar -rc libmuaompc.a *.o; echo "compiling $@"

clean:
	rm $(MPC_OBJ) $(MPC)libmuaompc.a mpcmydatmain_staticmem.o mpcmydatmain_staticmem
