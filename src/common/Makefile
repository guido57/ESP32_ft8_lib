# Native Linux build, including Raspberry Pi OS on a Raspberry Pi 5.
# The C decoder sources and C++ application sources must be compiled with
# their respective compilers, then linked with the C++ linker.

CC       ?= gcc
CXX      ?= g++
CPPFLAGS ?= -Isrc -DNATIVE_BUILD -DLOG_LEVEL=LOG_ERROR
CFLAGS   ?= -O3 -std=gnu11
CXXFLAGS ?= -O3 -std=gnu++17
LDLIBS   ?= -lm -latomic

TARGET := ft8_decoder
BUILD  := build/native

C_SOURCES := $(wildcard src/ft8/*.c src/fft/*.c) src/common/wave.c
CXX_SOURCES := src/main.cpp src/decoder_api.cpp src/decode_ft8.cpp src/subtract.cpp
OBJECTS := $(C_SOURCES:src/%.c=$(BUILD)/%.o) \
           $(CXX_SOURCES:src/%.cpp=$(BUILD)/%.o)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) -o $@ $^ $(LDLIBS)

$(BUILD)/%.o: src/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD) $(TARGET)
