CC = gcc
CFLAGS  = -O1 -Wall -Werror -lm -Wno-unused-result

TARGET = KalmanFilter

all: $(TARGET)

$(TARGET): $(TARGET).c
	$(CC) $(CFLAGS) -g -o $(TARGET) $(TARGET).c

clean:
	$(RM) $(TARGET)

