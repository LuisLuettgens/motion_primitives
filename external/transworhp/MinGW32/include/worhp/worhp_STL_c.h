#ifndef WORHP_STL_C_H
#define WORHP_STL_C_H

#include "C_std.h"

enum {
  QueueSize = 10
};

typedef struct WorhpQueueStruct {
  double val[QueueSize];
  int idx;
  int size;
} WorhpQueue;

DLL_PRIVATE void   InitWorhpQueue(WorhpQueue *const Q);
DLL_PRIVATE void   Push(WorhpQueue *const Q, double val);
DLL_PRIVATE double Peek(const WorhpQueue *const Q);
DLL_PRIVATE double Get (const WorhpQueue *const Q, int idx);

#endif
