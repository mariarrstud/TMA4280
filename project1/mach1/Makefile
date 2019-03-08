CC = mpicc
CFLAGS = -O2
TARGET = mach

$(TARGET): main.o mach.o vec_op.o
	$(CC) $(CFLAGS) -o $(TARGET) mach.o main.o vec_op.o -lm

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f mach mach.o main.o vec_op.o
