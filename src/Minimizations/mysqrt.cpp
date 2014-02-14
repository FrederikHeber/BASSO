#include "Minimizations.h"
#include "Table.h"

double mysqrt(const double a)
{
	double result = 0.;
	if (a >= 0) {
		if (((int)(a) == a) && (a < 10)) {
			result = sqrtTable[(int)(a)];
		} else {
#if defined (HAVE_LOG) && defined (HAVE_EXP)
			result = exp(log(x) * 0.5);
# else
			result = sqrt(a);
#endif
		}
	}
	return result;
}
