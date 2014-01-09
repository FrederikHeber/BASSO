#include "MathFunctions.h"

double mysqrt(const double a)
{
	double result = 0.;
	if (a >= 0)
		result = sqrt(a);
	return result;
}
