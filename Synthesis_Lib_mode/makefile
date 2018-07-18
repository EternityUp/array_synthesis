CC       :=  g++
CPPFLAGS :=  -g -std=c++11  -I /home/xzc/eigen3.3.4 -I ./inc/ -Wall -fPIC  -c 
AR       :=  ar
ARFLAGS  :=  crv
LD       :=  ld
LDFLAGS  :=  -shared -L ./libs/ -lArrayShapeAnalysis -lAudioFileProcess  -lFourierTransform -lBasicSignalProcess -lWavFileRW -lAgc -lm -lstdc++
RM       :=  rm -f 

SRC_PATH  = ./src/
INC_PATH  = ./inc/
DEMO_PATH = ./demo/
LIBS_DIR  = ./libs/

SRCS     := $(wildcard $(SRC_PATH)*.cpp)
OBJS     := $(patsubst %.cpp, %.o, $(SRCS))
LIB_SO   := $(LIBS_DIR)libArrayProcessSynthesis.so
LIB_A    := $(LIBS_DIR)libArrayProcessSynthesis.a
EXCS     := synthesis_lib_demo
LINKFLAGS := -g -std=c++11 -I /home/xzc/eigen3.3.4 -I ./inc/ -L $(LIBS_DIR) -Wl,-rpath=$(LIBS_DIR) -lArrayProcessSynthesis \
-lArrayShapeAnalysis -lAudioFileProcess -lWavFileRW -lAgc



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

$(OBJS):$(SRC_PATH)*.cpp
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
	
