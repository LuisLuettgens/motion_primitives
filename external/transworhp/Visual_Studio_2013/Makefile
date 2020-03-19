DIR = TransWORHP_Example
INC = $(DIR)/include

all: dist tar
#start:
#	cp ../example_release/transworhp.xml TransWORHP/release
#	cp ../example_release/worhp.xml TransWORHP/release


#install:
#	rm -r $(INC)/src
#	mkdir $(INC)/src
#	cp ../src/C_cs0.h ../src/TransWORHP.h ../src/Viewer.h ../src/butcher.h \
#	../src/diffstructure.h ../src/newdouble.h ../src/worhp_info.h \
##	../src/transworhpsharedmem.h $(INC)/src
#	rm -r $(INC)/gui
##	mkdir $(INC)/gui
#	mkdir $(INC)/gui/base
#	mkdir $(INC)/gui/glbase
#	mkdir $(INC)/gui/gui
#	mkdir $(INC)/gui/src_glplot
#	mkdir $(INC)/gui/toolbase
#	cp ../gui/base/color4.h $(INC)/gui/base
#	cp ../gui/base/vektor.h $(INC)/gui/base
#	cp ../gui/gui/sdlframe.h $(INC)/gui/gui
#	cp ../gui/gui/sdlscreen.h $(INC)/gui/gui
#	cp ../gui/gui/sdlcursor.h $(INC)/gui/gui
#	cp ../gui/gui/sdlthread.h $(INC)/gui/gui
#	cp ../gui/glbase/smoothmovement.h $(INC)/gui/glbase
#	cp ../gui/glbase/light.h $(INC)/gui/glbase
#	cp ../gui/glbase/font.h $(INC)/gui/glbase
#	cp ../gui/glbase/joystick.h $(INC)/gui/glbase
#	cp ../gui/glbase/ttfont.h $(INC)/gui/glbase
#	cp ../gui/glbase/texture.h $(INC)/gui/glbase
#	cp ../gui/glbase/viewport.h $(INC)/gui/glbase
#	cp ../gui/glbase/globject.h $(INC)/gui/glbase
#	cp ../gui/glbase/model.h $(INC)/gui/glbase
#	cp ../gui/base/point.h $(INC)/gui/base
#	cp ../gui/toolbase/toolmenu.h $(INC)/gui/toolbase
#	cp ../gui/toolbase/toolmenukey.h $(INC)/gui/toolbase
#	cp ../gui/toolbase/toolstatus.h $(INC)/gui/toolbase
#	cp ../gui/src_glplot/plot.h $(INC)/gui/src_glplot
#	cp ../gui/src_glplot/xopt_data.h $(INC)/gui/src_glplot
#	cp ../gui/src_glplot/functions.h $(INC)/gui/src_glplot
#	cp ../gui/src_xmlio/textout.h $(INC)
#	cp ../gui/src_xmlio/xmlio.h $(INC)
#	cp ../gui/src_xmlio/conversion.h $(INC)
#	rm -r $(DIR)/tutorial
#	mkdir $(DIR)/tutorial
#	cp ../tutorial/*.cpp $(DIR)/tutorial
#	rm -r $(DIR)/example
#	mkdir $(DIR)/example
#	cp ../example/*.cpp $(DIR)/example
#	cp ../example/*.c $(DIR)/example
#	cp TransWORHP/release/TransWORHP.lib $(DIR)/lib
#	cp TransWORHP/release/TransWORHP.dll $(DIR)/release
#	cp TransWORHP/release/transworhp.xml $(DIR)/release/
#	cp TransWORHP/release/worhp.xml $(DIR)/release/

#debug: $(DIR)/debug
#	cp TransWORHP/debug/TransWORHP.lib $(DIR)/lib
#	cp TransWORHP/debug/TransWORHP.dll $(DIR)/debug
#
#$(DIR)/debug:
#	mkdir $(DIR)/debug



dist:
	rm -f $(DIR)/release/*.pdb
	rm -f $(DIR)/release/*~
	rm -f $(DIR)/release/*.lib
	rm -f $(DIR)/release/*.exp
	rm -f $(DIR)/Visual?Studio/*.ncb
	rm -rf $(DIR)/Visual?Studio/*/Release
	rm -rf $(DIR)/Visual?Studio/*/Debug
	rm -f $(DIR)/Release/*.exp
	rm -f $(DIR)/Visual?Studio/Release/*.pdb
	rm -f $(DIR)/Visual?Studio/Release/*.lib
	rm -f $(DIR)/Visual?Studio/Release/*.exp




# Create archive file
DATE := $(shell date +"TransWORHP_Example_VC2010_%Y_%m_%d.tgz")

EXCLUDE = --exclude=*.vcproj.*.user --exclude=.svn

tar: 
	tar cvzf $(DATE) $(EXCLUDE) TransWORHP_Example
	cp $(DATE) /home/knauer/public_html/transworhp/download
	cp $(DATE) /home/knauer/public_html/amst

