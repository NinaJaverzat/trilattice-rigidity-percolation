#ifndef BASIC_H_
#define BASIC_H_

#include "defs.h"

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// INITALIZATION

void basic_init(Scalars* scalars, OrderParam* COP, OrderParam* ROP, OrderParam* PP);


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// GENERIC

void update_OP(OrderParam* OP, double val, double norm, int ind);               // Measure of step function

#endif /* BASIC_H_ */
