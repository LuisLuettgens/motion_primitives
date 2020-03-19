all:
	mkdir -p build
	cd build && cmake .. && make -j

init: get-worhp
	git submodule init
	git submodule update

init2: init
	mkdir -p build/example
	mkdir -p build/tutorial
	mkdir -p build/workshop
	mkdir -p build/Betts
	cp -n xml_files/*.xml build/example/
	cp -n xml_files/*.xml build/tutorial/
	cp -n xml_files/*.xml build/workshop/
	cp -n xml_files/*.xml build/Betts/

mini:
	mkdir -p build_mini
	cd build_mini && cmake -DBUILD-OPTIONS_TransWORHP_GRAPHICS=NO .. && make -j

# worhp mit "make libworhp" bauen
get-worhp:
	mkdir -p lib
	rm -f include/worhp/worhp/C_Worhp_Data.h
	cp -r ~/git/worhp/include/linux64/include/* include/worhp/
	cp -r ~/git/worhp/lib/linux64/* lib
	cp ~/git/worhp/bin/worhp.xml xml_files/

#get-worhp-titan: #muss angepasst werden
#	cp ../../../worhp/trunk/release/Linux/lib/include/worhp/* include/worhp/
#	cp -r ../../../worhp/trunk/release/Linux/lib/lib .

