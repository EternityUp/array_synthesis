GCC      :=  gcc
CFLAGS   :=  -g -std=c99 -Wall -fPIC  -c 
CC       :=  g++
CPPFLAGS :=  -g -std=c++11 -Wall -fPIC  -c 
AR       :=  ar
ARFLAGS  :=  crv
LD       :=  ld
LDFLAGS  :=  -shared -lm -lstdc++
RM       :=  rm -f 

SRC_PATH  = ./src/
INC_PATH  = ./inc/
DEMO_PATH = ./demo/
LIBS_DIR  = ./libs/

C_SRCS   := $(wildcard $(SRC_PATH)*.c)
CPP_SRCS := $(wildcard $(SRC_PATH)*.cpp)
C_OBJS   := $(patsubst %.c, %.o, $(C_SRCS))
CPP_OBJS := $(patsubst %.cpp, %.o, $(CPP_SRCS))
LIB_SO   := $(LIBS_DIR)libFourierTransform.so
LIB_A    := $(LIBS_DIR)libFourierTransform.a
EXCS     := fft_demo

LINKFLAGS := -std=c++11 -I ./inc/ -L $(LIBS_DIR) -Wl,-rpath=$(LIBS_DIR) -lFourierTransform -lm 



ifeq ($(LIB_TYPE), a)
	LINKFLAGS += -static
endif  

all:$(C_OBJS) $(CPP_OBJS) $(LIB_SO) $(LIB_A) $(EXCS)


$(EXCS):$(DEMO_PATH)*.cpp
	$(CC) $^ $(LINKFLAGS) -o $@ 
		

$(LIB_SO):$(C_OBJS) $(CPP_OBJS)
	$(LD) $(LDFLAGS) $^ -o $@

$(LIB_A):$(C_OBJS) $(CPP_OBJS)
	$(AR) $(ARFLAGS) $@ $^

$(C_OBJS):$(SRC_PATH)*.c 
	$(GCC) $(CFLAGS) $^ -I $(INC_PATH) -o $@
	
$(CPP_OBJS):$(SRC_PATH)*.cpp
	$(CC) $(CPPFLAGS) $^ -I $(INC_PATH) -o $@

	
.PHONY:clean
clean:
	$(RM) $(LIB_SO) $(LIB_A) $(C_OBJS) $(CPP_OBJS) $(EXCS)

	
.PHONY:print
print:
	@echo CFLAGS    = $(CFLAGS)
	@echo CPPFLAGS    = $(CPPFLAGS)
	@echo ARFLAGS     = $(ARFLAGS)
	@echo LDFLAGS     = $(LDFLAGS)
	@echo LINKFLAGS   = $(LINKFLAGS)
	@echo C_OBJS      = $(C_OBJS)
	@echo CPP_OBJS    = $(CPP_OBJS)
	@echo C_SRCS      = $(C_SRCS)
	@echo CPP_SRCS    = $(CPP_SRCS)
	@echo LIB_TYPE    = $(LIB_TYPE)
	
