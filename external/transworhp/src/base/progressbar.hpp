#include <time.h>

class ProgressBar {
public:
	ProgressBar(unsigned int maxval_);
	~ProgressBar();
	void status(unsigned int current);

private:
	void set_buffer(float percent);
	float time_elapsed();
	clock_t starttime;
	char buffer[80];
	unsigned int size;
	unsigned int maxval;
};
