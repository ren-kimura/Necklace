CC = gcc
CFLAGS = -Wall -Wextra -O2 -g
LDFLAGS = -lm

TARGETS = nkl vfa

BASE_OBJS = ds.o utils.o cov.o stat.o write.o veri.o out.o eutils.o

NKL_OBJS = nkl.o $(BASE_OBJS)
VFA_OBJS = vfa.o ds.o utils.o stat.o veri.o eutils.o out.o

all: $(TARGETS)

nkl: $(NKL_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

vfa: $(VFA_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o $(TARGETS)

.PHONY: all clean