#pragma once

/** Erzwinge Kompatibilitaet unter Win zwischen Debug- und Release-Versionen */
//#define _ITERATOR_DEBUG_LEVEL 0
//#define _HAS_ITERATOR_DEBUGGING 0

/** Art der Diskretisierung */
enum class TWdiscretizationType {Trapez=0, HermiteSimpson=1, Euler=2, Lobatto=3, MultipleShooting=10, Pseudospectral=11};
