CC       :=  g++
CPPFLAGS :=  -g -std=c++11 -Wall -fPIC  -c 
RM       :=  rm -f 

SRC_PATH   =  ./src/
INC_PATH1  =  ./inc/
INC_PATH2  =  /home/xzc/eigen3.3.4/
DEMO_PATH  =  ./demo/
LIBS_DIR   =  ./libs/

SRCS     := $(wildcard $(SRC_PATH)*.cpp)
OBJS     := $(patsubst %.cpp, %.o, $(SRCS))
EXCS     := synthesis_demo

LINKFLAGS := -g -Wall -fPIC -std=c++11 -I $(INC_PATH1) -I $(INC_PATH2) -L $(LIBS_DIR) -lArrayShapeAnalysis \
-lFourierTransform -lAudioFileProcess -lBasicSignalProcess -lWavFileRW



all:$(OBJS) $(EXCS)


$(EXCS):$(DEMO_PATH)*.cpp
	$(CC) $^ $(LINKFLAGS) -o $@	

%.o:%.cpp
	$(CC) $(CPPFLAGS) $^ -I $(INC_PATH1) -I $(INC_PATH2) -o $@	

.PHONY:clean
clean:
	$(RM) $(OBJS) $(EXCS)

.PHONY:print	
print:
	@echo CPPFLAGS  = $(CPPFLAGS)
	@echo LINKFLAGS = $(LINKFLAGS)
	@echo OBJS      = $(OBJS)
	@echo SRCS      = $(SRCS)
	
