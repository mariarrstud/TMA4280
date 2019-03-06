CC = gcc
CFLAGS = -O2
TARGET = mach
LIB = libutils.a
MAINSRC = main.c

SRC = $(wildcard *.c)

OBJ = $(patsubst %.c, %.o, $(SRC))

MAINOBJ = $(patsubst %.c, %.o, $(MAINSRC))
LIBOBJ = $(patsubst $(MAINOBJ),,$(OBJ))

%.o:%.c
	$(CC) $(CFLAGS) -c $<

all: $(LIB) $(MAINOBJ)
	$(CC) -o $(TARGET) $(MAINOBJ) $(LIB) -lm

$(LIB): $(LIBOBJ)
	@ar r $@ $(LIBOBJ)

clean:
	rm -f $(TARGET) $(LIB) $(OBJ)

list_libobjs: $(LIB)
	@echo "All the objects contained in the static library are printed below:\n"
	nm $(LIB)
