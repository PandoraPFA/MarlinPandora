#Path to project directory
PROJECT_DIR = YOUR_PATH_HERE

#Paths to project dependencies, note monitoring is optional
GEAR_DIR = YOUR_PATH_HERE
LCIO_DIR = YOUR_PATH_HERE
MARLIN_DIR = YOUR_PATH_HERE
MARLINUTIL_DIR = YOUR_PATH_HERE
PANDORAPFANEW_DIR = YOUR_PATH_HERE
PANDORAMONITORING_DIR = YOUR_PATH_HERE

DEFINES = -DUSE_GEAR=1
ifdef MONITORING
    DEFINES  += -DMONITORING=1
endif

PROJECT_INCLUDE_DIR = $(PROJECT_DIR)/include/
PROJECT_SOURCE_DIR  = $(PROJECT_DIR)/src/
PROJECT_LIBRARY_DIR = $(PROJECT_DIR)/lib/

INCLUDES  = -I$(PROJECT_INCLUDE_DIR)
INCLUDES += -I$(GEAR_DIR)/include/
INCLUDES += -I$(LCIO_DIR)/include/
INCLUDES += -I$(MARLIN_DIR)/include/
INCLUDES += -I$(MARLINUTIL_DIR)/include/
INCLUDES += -I$(PANDORAPFANEW_DIR)/Framework/include/
INCLUDES += -I$(PANDORAPFANEW_DIR)/FineGranularityContent/include/ -I$(PANDORAPFANEW_DIR)/KMeansContent/include/
ifdef MONITORING
    INCLUDES += -I$(PANDORAMONITORING_DIR)/include/
endif

CC = g++
CFLAGS = -c -Wall -g -w -fPIC -O2
CFLAGS += $(INCLUDES)
ifdef BUILD_32BIT_COMPATIBLE
    CFLAGS += -m32
endif

SOURCES = $(wildcard $(PROJECT_SOURCE_DIR)*.cc)

OBJECTS = $(SOURCES:.cc=.o)

LIBS = -L$(GEAR_DIR)/lib -lgear
LIBS += -L$(LCIO_DIR)/lib -llcio
LIBS += -L$(MARLIN_DIR)/lib -lMarlin
LIBS += -L$(MARLINUTIL_DIR)/lib -lMarlinUtil
LIBS += -L$(PANDORAPFANEW_DIR)/lib -lPandoraFramework -lPandoraFineGranularityContent -lPandoraKMeansContent
ifdef BUILD_32BIT_COMPATIBLE
    LIBS += -m32
endif

LDFLAGS = $(LIBS) -Wl,-rpath

LIBRARY = $(PROJECT_LIBRARY_DIR)/libMarlinPandora.so

all: $(SOURCES) $(OBJECTS)
	$(CC) $(OBJECTS) $(LIBS) -shared -o $(LIBRARY)

$(LIBRARY): $(OBJECTS)
	$(CC) $(LDFLAGS) -fPIC $(OBJECTS) -o $@

.cc.o:
	$(CC) $(CFLAGS) $(DEFINES) $< -o $@

clean:
	rm -f $(OBJECTS)
	rm -f $(LIBRARY)
