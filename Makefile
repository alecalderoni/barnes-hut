CC := gcc-15
CFLAGS := -Wall -O3 -Iinclude -fopenmp
LDFLAGS := -lm -fopenmp

SRC := src/main.c src/galactic_dynamics.c
OBJ := $(SRC:.c=.o)
BIN := app

.PHONY: all clean

all: $(BIN)

$(BIN): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDFLAGS)

src/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

run:
	caffeinate ./app

clean:
	rm -f src/*.o $(BIN) dati/*.dat
	
