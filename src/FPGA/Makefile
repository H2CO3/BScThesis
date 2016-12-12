SYNTHESIS ?= 0
NDEBUG ?= 0

OPSYS = $(shell uname | tr '[[:upper:]]' '[[:lower:]]')

ifeq ($(OPSYS), linux)
	CXX = g++
else
	CXX = xcrun -sdk macosx clang++
endif

LD = $(CXX)

CXFLAGS = -std=c++11 -c -Iinclude -O3 -flto \
	-Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-label

LDFLAGS = -O3 -flto

ifneq ($(SYNTHESIS), 0)
	CXFLAGS += -D__SYNTHESIS__
else
	CXFLAGS += -U__SYNTHESIS__
endif

ifneq ($(NDEBUG), 0)
	CXFLAGS += -DNDEBUG
else
	CXFLAGS += -UNDEBUG
endif

all: clean align

align: align.o main.o
	$(LD) $(LDFLAGS) -o $@ $^

%.o:%.cc
	$(CXX) $(CXFLAGS) -o $@ $<

clean:
	rm -f align align *.o

.PHONY: all clean
