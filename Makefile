TARGET_EXEC ?= worker

SRCS := $(shell find . -name '*.c')
OBJECTS := $(SRCS:%.c=%.o)

# Add include dir and lib path, with all the warnings
CFLAGS=-c -Wall -Werror -O2
LDFLAGS=-lgmp

$(TARGET_EXEC): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

# c source
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJECTS)
	$(RM) $(TARGET_EXEC)