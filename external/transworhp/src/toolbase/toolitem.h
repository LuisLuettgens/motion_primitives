#pragma once

#include "toolrect.h"

namespace tw {


typedef void (*signal1) ();

/** @ingroup toolbase
 *  @brief Item
 */
class ToolItem {
public:
	ToolItem();
	ToolItem(const ToolRect &rect_);
	virtual ~ToolItem();

	bool contains(int x_,int y_);

	virtual void Draw(bool isActive, const Point<int> &mouse);

	virtual void MouseClick(int button, int state, const Point<int> &p);
	virtual void Key(SDL_Keysym &Keysym);
	//virtual void SpecialKey(Uint8 *keys);
	virtual bool TakesKey();
	virtual void Apply();
	virtual void Discard();
	void Signal(signal1 s_);

//protected:
	ToolRect rect;
	signal1 s;
};


/** @ingroup toolbase
 *  @brief Label
 */
template <class T>
class ToolLabel : public ToolItem {
public:
//	ToolLabel(T &r) : ref(r) {};
	ToolLabel(const ToolRect &rect_, T &r);

	void Draw(bool isActive, const Point<int> &mouse);
	void Update();
	// void Set(T &r);

	T &ref;
};

/** @ingroup toolbase
 *  @brief Label
 */
template <class T>
class ToolFrame : public ToolItem {
public:
//	ToolLabel(T &r) : ref(r) {};
	ToolFrame(const ToolRect &rect_, T &r);

	void Draw(bool isActive, const Point<int> &mouse);
	void Update();
	// void Set(T &r);

	T &ref;
	ToolRect r2;

};

/** @ingroup toolbase
 *  @brief Template fuer verschiedene Speichertypen: Werte uebernehmen oder vergessen
 */
template <class T>
class ToolData : public ToolItem {
public:
	ToolData(const ToolRect &rect_, T &r)
			:  ToolItem(rect_),ref(r), value(r) {}

	virtual ~ToolData() {}

	virtual void Apply() {
		ref = value;
		//std::cout << "Item: " << ToString(ref) << std::endl;
	}
	virtual void Discard() {

		//std::cout << ToString(ref) <<" to " << ToString(value) << std::endl;
		value = ref;
	}

	T &ref;
	T value;
};

}
