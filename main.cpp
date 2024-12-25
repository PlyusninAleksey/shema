#include <iostream>
#include <cmath>
#include "scheme.h"
#include <vector>
using std::vector;

int main()
{
	int n = 1000;
	Scheme shema;
	shema.calculate_test(n);
	for (int i = 0; i < n+1; i++)
	{
		std::cout << shema.x[i] << "\t" << shema.v[i] << "\t" << shema.v2[i] << "\t" << shema.diff[i] << "\n";
	}
	return 0;
}