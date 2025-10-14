CC = gcc
CFLAGS = -Wall -Wextra -O2 -g
LDFLAGS = -lm

TARGETS = nkl eu

BASE_OBJS = ds.o utils.o cov.o stat.o write.o veri.o

NKL_OBJS = nkl.o $(BASE_OBJS) out.o
EU_OBJS = eu.o $(BASE_OBJS) eutils.o

all: $(TARGETS)

nkl: $(NKL_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

eu: $(EU_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o $(TARGETS)

.PHONY: all clean