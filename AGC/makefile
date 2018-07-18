CC       :=  g++
CPPFLAGS :=  -g -std=c++11 -Wall -fPIC  -c 
AR       :=  ar
ARFLAGS  :=  crv
LD       :=  ld
LDFLAGS  :=  -shared -lm
RM       :=  rm -f 

SRC_PATH  = ./src/
INC_PATH  = ./inc/
DEMO_PATH = ./demo/
LIBS_DIR  = ./libs/

SRCS     := $(wildcard $(SRC_PATH)*.c)
OBJS     := $(patsubst %.c, %.o, $(SRCS))
LIB_SO   := $(LIBS_DIR)libAgc.so
LIB_A    := $(LIBS_DIR)libAgc.a
EXCS     := agc_demo

#LINKFLAGS := -D__SIMULATE_POSITIONS__ -I ./inc/ -L. -Wl,-rpath=./ -lArrayShapeAnalysis -lm 
LINKFLAGS := -I ./inc/ -L $(LIBS_DIR) -Wl,-rpath=$(LIBS_DIR) -lAgc -lm 



ifeq ($(LIB_TYPE), a)
	LINKFLAGS += -static
endif  

all:$(OBJS) $(LIB_SO) $(LIB_A) $(EXCS)


$(EXCS):$(DEMO_PATH)*.c
	$(CC) $^ $(LINKFLAGS) -o $@ 
		

$(LIB_SO):$(OBJS)
	$(LD) $(LDFLAGS) $^ -o $@

$(LIB_A):$(OBJS)
	$(AR) $(ARFLAGS) $@ $^

%.o:%.c
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
	
