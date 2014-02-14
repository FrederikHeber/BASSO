// A simple program that computes the square root of a number
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "BassoConfig.h"
#ifdef USE_MYMATH
#include "Minimizations/Minimizations.h"
#endif

int main (int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stdout,"%s version %d.%d\n",
			argv[0],
			Basso_VERSION_MAJOR,
			Basso_VERSION_MINOR);
		fprintf(stdout, "Usage: %s number\n",
			argv[0]);
		return 1;
	}
	double inputValue = atof(argv[1]);

#ifdef USE_MYMATH
	double outputValue = mysqrt(inputValue);
#else
	double outputValue = sqrt(inputValue);
#endif
	if (inputValue < 0)
		outputValue = 0.;

	fprintf(stdout,"The square root of %g is %g\n",
		inputValue, 
		outputValue);
	return 0;
}
