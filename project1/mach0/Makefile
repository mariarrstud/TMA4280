CC = gcc
CFLAGS = -O2
TARGET = mach

$(TARGET): main.o mach.o
	$(CC) $(CFLAGS) -o $(TARGET) mach.o main.o -lm

utest: utest.o mach.o
	$(CC) $(CFLAGS) -o utest utest.o mach.o -lm

vtest: vtest.o mach.o
	$(CC) $(CFLAGS) -o vtest vtest.o mach.o -lm

test: utest vtest
	./utest ./vtest

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f utest utest.o vtest vtest.o mach.o mach main.o
