CC = gcc
CFLAGS = -O -Wall
LDFLAGS = -lm

TARGET = tmzdemo2.exe

SRC = tmzdemo2.c \
      fdtd/abctmz.c \
      fdtd/grid1dez.c \
      fdtd/gridtmz.c \
      fdtd/ricker.c \
      fdtd/solver.c \
      fdtd/solver_tr.c \
      fdtd/tfsftmz.c \
      fdtd/pml.c \
      fdtd/pml_tr.c \
      fdtd/updatetmz.c \
      fdtd/updatetmz_tr.c \
      utility/fdtd-config.c \
      utility/fdtd-io.c \
      utility/snapshot2d.c

OBJ = $(SRC:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(OBJ) $(LDFLAGS) -o $@
	powershell -NoProfile -Command "Get-ChildItem -Recurse -Filter '*.o' -File | Remove-Item -Force"

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	powershell -NoProfile -Command "Get-ChildItem -Recurse -Filter '*.o' -File | Remove-Item -Force; Remove-Item -Force -ErrorAction SilentlyContinue '$(TARGET)','ez_field.csv'; Remove-Item -Force -ErrorAction SilentlyContinue sim.*"

rebuild: clean all

.PHONY: all clean rebuild
