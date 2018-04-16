# for C++ define  CC = g++
CC = g++
CFLAGS = -g -Wall -fPIC -m64 -std=gnu++14
LFLAGS = -L$(ROOTSYS)/lib 
INC =	-I$(ROOTSYS)/include -I$(O2_ROOT)/include -I$(FAIRROOT_ROOT)/include -I./
TGT =	libMatBud.so
DICT=	MatBudDict.cxx
DICTO=	MatBudDict.o

SRC = 	Ray.cxx MatLayerCyl.cxx MatLayerCylSet.cxx

HDR =	$(SRC:.cxx=.h) 

OBJ = 	$(SRC:.cxx=.o)


.PHONY: depend clean

all: 	$(TGT)
	@echo creating libMatBud.so

$(TGT):	$(OBJ) $(DICTO)
	$(CC) $(CFLAGS)  -shared -o $(TGT) $(OBJ) $(DICTO) `root-config --ldflags` $(LFLAGS)

# pull in dependency info for *existing* .o files
-include $(OBJ:.o=.d)

%.o : %.cxx
	$(CC) $(CFLAGS) $(INC) -c $<  -o $@
	$(CC) -MM $(CFLAGS) $(INC) -c $*.cxx > $*.d
	@cp -f $*.d $*.d.tmp
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp

clean:
	rm -f *.o *~ *.so *.d *Dict.{h,cxx}

$(DICT): $(HDR) MatBudLinkDef.h
	rootcint -f $@ -c $(INC) $(HDR) $^


depend: $(SRC)
	makedepend $(INC) $^

# DO NOT DELETE THIS LINE -- make depend needs it
