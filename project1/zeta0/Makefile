CC = gcc
CFLAGS = -O2
TARGET = zeta

$(TARGET): main.o zeta.o
	$(CC) $(CFLAGS) -o $(TARGET) zeta.o main.o -lm

utest: utest.o zeta.o
	$(CC) $(CFLAGS) -o utest utest.o zeta.o -lm

vtest: vtest.o zeta.o
	$(CC) $(CFLAGS) -o vtest vtest.o zeta.o -lm

test: utest vtest
	./utest ./vtest

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f utest utest.o vtest vtest.o zeta zeta.o main.o
