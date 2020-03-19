DIR = TransWORHP
INC = $(DIR)/include
LIB = $(DIR)/lib
WORHP = $(INC)/worhp
WORHPL = $(INC)/worhp/lib


install: 
	mkdir -p $(DIR)
	mkdir -p $(INC)
	mkdir -p $(LIB)
	mkdir -p $(WORHP)

	mkdir -p $(DIR)/example_release
	mkdir -p $(DIR)/example_release/modell
	mkdir -p $(DIR)/tutorial_release

	cp mainMakefile $(DIR)/Makefile

	rm -fr $(INC)/core
	mkdir $(INC)/core
	cp ../src/core/TransWORHP.h $(INC)/core
	cp ../src/core/Viewer.h $(INC)/core
	cp ../src/core/butcher.h $(INC)/core
	cp ../src/core/timing.h $(INC)/core
	cp ../src/core/MagicTransWORHP.h $(INC)/core
	cp ../src/core/ExplTransWORHP.h $(INC)/core
	cp ../src/core/PmTransWORHP.h $(INC)/core
	cp ../src/core/LobattoPmTransWORHP.h $(INC)/core
	cp ../src/core/FullDisTransWORHP.h $(INC)/core
	cp ../src/core/diffstructure.h $(INC)/core
	cp ../src/core/newdouble.h $(INC)/core
	cp ../src/core/worhp_info.h $(INC)/core
	cp ../src/core/TWsharedmemory.h $(INC)/core
	cp ../src/core/TWparameter.h $(INC)/core
	cp ../src/core/TWcount.h $(INC)/core
	cp ../src/core/TWfolder.h $(INC)/core
	cp ../src/core/TWconsole.h $(INC)/core
	cp ../src/core/TWdebug.h $(INC)/core
	cp ../src/core/TWproblem.h $(INC)/core
	cp ../src/core/twstatus.h $(INC)/core
	cp ../src/core/TWGUIconfig.h $(INC)/core

	rm -rf $(INC)/base
	mkdir $(INC)/base
	cp ../src/base/defines.h $(INC)/base
	cp ../src/base/color4.h $(INC)/base
	cp ../src/base/vektor.h $(INC)/base
	cp ../src/base/point.h $(INC)/base

	rm -rf $(INC)/gui
	mkdir $(INC)/gui
	cp ../src/gui/sdlframe.h $(INC)/gui
	cp ../src/gui/sdlscreen.h $(INC)/gui
	cp ../src/gui/sdlcursor.h $(INC)/gui
	cp ../src/gui/sdlthread.h $(INC)/gui

	rm -rf $(INC)/glbase
	mkdir $(INC)/glbase
	cp ../src/glbase/smoothmovement.h $(INC)/glbase
	cp ../src/glbase/light.h $(INC)/glbase
	cp ../src/glbase/font.h $(INC)/glbase
	cp ../src/glbase/joystick.h $(INC)/glbase
	cp ../src/glbase/texture.h $(INC)/glbase
	cp ../src/glbase/viewport.h $(INC)/glbase
	cp ../src/glbase/globject.h $(INC)/glbase
	cp ../src/glbase/model.h $(INC)/glbase

	rm -rf $(INC)/toolbase
	mkdir $(INC)/toolbase
	cp ../src/toolbase/toolmenu.h $(INC)/toolbase
	cp ../src/toolbase/toolmenukey.h $(INC)/toolbase
	cp ../src/toolbase/toolstatus.h $(INC)/toolbase

	rm -rf $(INC)/glplot
	mkdir $(INC)/glplot
	cp ../src/glplot/plot.h $(INC)/glplot
	cp ../src/glplot/baseplot.h $(INC)/glplot
	cp ../src/glplot/threedplot.h $(INC)/glplot
	cp ../src/glplot/xopt_data.h $(INC)/glplot
	cp ../src/glplot/functions.h $(INC)/glplot
	
	rm -rf $(INC)/xmlio
	mkdir $(INC)/xmlio
	cp ../src/xmlio/src/xmlio/textout.h $(INC)/xmlio
	cp ../src/xmlio/src/xmlio/xmlio.h $(INC)/xmlio
	cp ../src/xmlio/src/xmlio/conversion.h $(INC)/xmlio

	rm -fr $(DIR)/tutorial
	mkdir $(DIR)/tutorial
	cp ../tutorial/*.cpp $(DIR)/tutorial
	cp tutorialcmakelists.txt $(DIR)/tutorial/CMakeLists.txt
	
#	rm -fr $(DIR)/tutorial-mini
#	mkdir $(DIR)/tutorial-mini
#	cp tutorialminicmakelists.txt $(DIR)/tutorial-mini/CMakeLists.txt

	rm -rf $(DIR)/example
	mkdir $(DIR)/example
	cp ../example/*.cpp $(DIR)/example
	cp ../example/*.c $(DIR)/example
	cp ../example/*.h $(DIR)/example
	cp examplecmakelists.txt $(DIR)/example/CMakeLists.txt

	cp -rf ../include/worhp/worhp/* $(WORHP)

	cp ../build/src/*.so $(LIB)
	cp ../build_mini/src/*.so $(LIB)

	cp ../lib/*worhp.so* $(LIB)

	cp ../xml_files/*.xml $(DIR)/example_release
	cp ../example_release/modell/* $(DIR)/example_release/modell
	cp ../xml_files/*.xml $(DIR)/tutorial_release




# Create archive file
DATE := $(shell date +"TransWORHP_Linux64_%Y_%m_%d.tgz")

EXCLUDE = --exclude=build --exclude=.svn

tar:
	tar cvzf $(DATE) $(EXCLUDE) TransWORHP
#	cp $(DATE) /home/knauer/public_html/transworhp/download


ftp:
	ftpscript $(DATE)

