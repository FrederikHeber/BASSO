#include "MathFunctions.h"

double mysqrt(const double a)
{
	double result = 0.;
	if (a >= 0) {
#if defined (HAVE_LOG) && defined (HAVE_EXP)
		result = exp(log(x) * 0.5);
# else
		result = sqrt(a);
#endif
	}
	return result;
}
