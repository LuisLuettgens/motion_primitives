include ../Makefile.inc
include ../../Makefile.inc


CPPFILES := $(wildcard *.cpp)
APPS := $(subst .cpp,,$(CPPFILES))

all: $(APPS)

DYNLIBS = -L$(LIBDIR) -l$(XMLIO)

%: %.cpp $(LIBS)
	$(GCC) $(CFLAGS) -I$(INCLUDEDIR) $< $(DYNLIBS) -o $@

test: xmlview
	xmlview test.xml ASCII
	xmlview test.xml CONSOLE
	xmlview test.xml RTF > test.rtf
	xmlview test.xml HTML > test.html
	xmlview test.xml CPP
	xmlview /home/knauer/include/font/font.xml CPP > xmlfont.h


#	xmltext CPP

clean:
	rm -f *.cpp~ *.h~ $(APPS) *.html *.rtf

