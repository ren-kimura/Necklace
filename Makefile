CC = gcc
CFLAGS = -Wall -Wextra -O2 -g
LDFLAGS = -lm

SRCS = necklace.c ds.c utils.c mbm.c rep_utils.c
OBJS = $(SRCS:.c=.o)

TARGET = necklace

all: $(TARGET)

$(TARGET): $(OBJS)
		$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

%.o: %.c
		$(CC) $(CFLAGS) -c $< -o $@

clean:
		rm -f $(OBJS) $(TARGET)

.PHONY: all clean