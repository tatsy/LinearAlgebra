#include <ctime>
#include "MTRand.h"

extern "C" {
#include "mt19937ar.h"
}

MTRand::MTRand()
{
	init_genrand((unsigned long)time(0));
}

int MTRand::randInt(const int n) const
{
	return genrand_int31() % n;
}

double MTRand::randReal() const
{
	return genrand_real2();
}
