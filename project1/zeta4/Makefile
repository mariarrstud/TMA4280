CC = mpicc
CFLAGS = -fopenmp
TARGET = zeta

$(TARGET): main.o zeta.o vec_op.o
	$(CC) $(CFLAGS) -o $(TARGET) zeta.o main.o vec_op.o -lm

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f zeta zeta.o main.o vec_op.o
