CC       :=  g++
CPPFLAGS :=  -g -std=c++11 -DWEBRTC_ARCH_LITTLE_ENDIAN -Wall -fPIC -c 
#CPPFLAGS :=  -g -std=c++11 -DWEBRTC_ARCH_LITTLE_ENDIAN -Wall -fPIC -fno-use-cxa-atexit -c 
AR       :=  ar
ARFLAGS  :=  crv
LD       :=  ld
LDFLAGS  :=  -shared
RM       :=  rm -f 

SRC_PATH  = ./src/
INC_PATH  = ./inc/
DEMO_PATH = ./demo/
LIBS_DIR  = ./libs/

SRCS     := $(wildcard $(SRC_PATH)*.cpp)
OBJS     := $(patsubst %.cpp, %.o, $(SRCS))
LIB_SO   := $(LIBS_DIR)libWavFileRW.so
LIB_A    := $(LIBS_DIR)libWavFileRW.a
EXCS     := basic_signal_process_lib_demo

LINKFLAGS := -g -std=c++11 -I ./inc/ -L $(LIBS_DIR) -Wl,-rpath=$(LIBS_DIR) -lWavFileRW -lm 



ifeq ($(LIB_TYPE), a)
	LINKFLAGS += -static
endif  

all:$(OBJS) $(LIB_SO) $(LIB_A) $(EXCS)


$(EXCS):$(DEMO_PATH)*.cpp
	$(CC) $^ $(LINKFLAGS) -o $@ 
		

$(LIB_SO):$(OBJS)
	$(LD) $(LDFLAGS) $^ -o $@

$(LIB_A):$(OBJS)
	$(AR) $(ARFLAGS) $@ $^

%.o:%.cpp
	$(CC) $(CPPFLAGS) $^ -I $(INC_PATH) -o $@	

.PHONY:clean
clean:
	$(RM) $(LIB_SO) $(LIB_A) $(OBJS) $(EXCS)

.PHONY:print	
print:
	@echo CPPFLAGS  = $(CPPFLAGS)
	@echo ARFLAGS   = $(ARFLAGS)
	@echo LDFLAGS   = $(LDFLAGS)
	@echo LINKFLAGS = $(LINKFLAGS)
	@echo OBJS      = $(OBJS)
	@echo SRCS      = $(SRCS)
	@echo LIB_TYPE  = $(LIB_TYPE)
	
