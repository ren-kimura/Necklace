CC = gcc
CFLAGS = -Wall -Wextra -O2 -g
LDFLAGS = -lm

SRCS = nkl.c ds.c utils.c cov.c out.c stat.c
OBJS = $(SRCS:.c=.o)

TARGET = nkl

all: $(TARGET)

$(TARGET): $(OBJS)
		$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

%.o: %.c
		$(CC) $(CFLAGS) -c $< -o $@

clean:
		rm -f $(OBJS) $(TARGET)

.PHONY: all clean