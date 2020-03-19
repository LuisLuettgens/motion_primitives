#pragma once

#include "worhp/worhp.h"

#include <iostream>

std::ostream &operator<<(std::ostream& os, const OptVar &o);

void print(WorhpMatrix &M, int index, std::ostream &os = std::cout);
