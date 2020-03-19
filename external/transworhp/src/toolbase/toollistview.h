//
// C++ Interface: toolbox
//
// Description:
//
//
// Author: Matthias Knauer,,, <tulio@visurgis>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TOOLLISTVIEW_H
#define TOOLLISTVIEW_H

#include "tool.h"

#include <vector>

#include "toolwindow.h"
#include "toolrect.h"


/**
 @author Matthias Knauer,,, <tulio@visurgis>
*/

class ToolListitem {
public:
	ToolListitem() {}
	virtual ~ToolListitem() {}
	virtual const std::string GetString() const = 0;
	virtual void Draw(ToolRect &r) {}
	;
};

/** Show data, no update possible */
template <class T>
class ToolListD : public ToolListitem {
public:
	ToolListD(T v, int m=0) : value(v), mode(m) {}
	T value;
	const std::string GetString() const;
	void Draw(ToolRect &r);
	int mode;
};

/** Show reference to data, update automatic */
template <class T>
class ToolListR : public ToolListitem {
public:
	ToolListR(T& v, int m=0) : ref(v), mode(m) {}
	T &ref;
	const std::string GetString() const;
	void Draw(ToolRect &r);
	int mode;
};

template <class T>
class ToolListC : public ToolListD<T> {
public:
	ToolListC(T v, color4 c) : ToolListD<T>(v), col(c) {}

	void Draw(ToolRect &r);
	color4 col;
};

class ToolListdata {
public:
	ToolListdata(int l=0) : open(0), level(l) {}
	~ToolListdata() {
		clearpointer(data);
	}

	void AddItem(ToolListitem*);

	std::vector<ToolListitem *> data;
	bool open;
	int level;
};


class ToolListcolumn {
public:
	ToolListcolumn(const std::string &s, int w) : text(s), width(w) {}
	~ToolListcolumn() {}

	//void AddItem(ToolListitem*);

	//std::vector<ToolListitem *> data;

	std::string text;
	int width;
};


class ToolListview : public ToolWindow {
public:
	ToolListview(const std::string &t, int w, int h, bool write=true);
	virtual ~ToolListview();

	void Draw(bool isActive, int HT, const Point<int> &mouse, bool border=true);
	void Update() {}
	void Init() {}

	void WriteText();
	void Add(ToolListdata *a) {
		list.push_back(a);
	}
	void AddColumn(const std::string &s, int w);
	//protected:
	bool mousebutton(int mx, int my, int button, int type);

	void SetMessage(int id) {
		message_id=id;
	}
protected:
	std::vector<ToolListdata *> list;
	std::vector<ToolListcolumn *> titles;

	int message_id;
};


#endif
