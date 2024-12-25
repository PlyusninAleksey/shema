#pragma once
#include <cmath>
#include <vector>
using std::vector;

const double M_PI = 3.14159265358979323846;

double k_main(double x, double xi)
{
	if (x <= xi)
	{
		return sqrt(2) * cos(x);
	}
	else
	{
		return 2;
	}
}

double q_main(double x, double xi)
{
	if (x <= xi)
	{
		return x;
	}
	else
	{
		return x * x;
	}
}

double f_main(double x, double xi)
{
	if (x <= xi)
	{
		return sin(2.0 * x);
	}
	else
	{
		return sin(x);
	}
}

double k_test(double x, double xi)
{
	if (x <= xi)
	{
		return 1.0;
	}
	else
	{
		return 2.0;
	}
}

double q_test(double x, double xi)
{
	if (x <= xi)
	{
		return M_PI / 4.0;
	}
	else
	{
		return pow((M_PI / 4.0), 2);
	}
}

double f_test(double x, double xi)
{
	if (x <= xi)
	{
		return 1.0;
	}
	else
	{
		return sqrt(2.0) / 2;
	}
}



vector<double> solveMatrix(vector<double> A, vector<double>& C, vector<double> B, vector<double> Fi, double mu1, double mu2, int n)
{
	vector<double> y(n + 1);
	vector<double> alpha(n);
	vector<double> beta(n);

	alpha[0] = 0.0;
	beta[0] = mu1;
	for(int i = 0; i < n-1; i++)
	{
		alpha[i + 1] = B[i] / (C[i] - alpha[i] * A[i]);
		beta[i + 1] = (Fi[i] + A[i] * beta[i]) / (C[i] - alpha[i] * A[i]);
	}

	y[n] = mu2;

	for(int i = n - 1; i >= 0; i--)
	{
		y[i] = alpha[i] * y[i + 1] + beta[i];
	}
	return y;
}

vector<double> calc_true_sol(double levGran, double pravGran, double xi, double n)
{
	const double C1 = -0.2512338358800897;
	const double C2 = -0.5493987964341227;
	const double C3 = -0.09925199591055987;
	const double C4 = 0.046413492642686384;

	vector<double> trueSol(n+1);
	double h = (pravGran - levGran) / n;
	int i = 0;
	 

	while(levGran + h*i < xi)
	{
		trueSol[i] = C1 * exp((sqrt(M_PI) / 2.0) * (levGran + h * i)) + C2 * exp((-sqrt(M_PI) / 2.0) * (levGran + h * i)) + (4.0 * sqrt(2.0)) / M_PI;
		i++;
	}

	for (; i < n + 1; i++)
	{
		trueSol[i] = C3 * exp((M_PI / (4.0 * sqrt(2.0))) * (levGran + h * i)) + C4 * exp((-M_PI / (4.0 * sqrt(2.0))) * (levGran + h * i)) + (8.0 * sqrt(2.0)) / (M_PI * M_PI);
	}
	return trueSol;

}

class Scheme
{
public:
	double mu1;
	double mu2;
	double xi; //кси
	double levGran;
	double pravGran;
	vector<double> x; //иксы
	vector<double> v; //численное решение
	vector<double> v2; //либо истинное, либо основное с большой сеткой
	vector<double> diff; //разница между vi и v2i
	Scheme(double mu1 = 1.0, double mu2 = 1.0, double xi = M_PI / 4, double a = 0.0, double b = 1.0) :mu1(mu1), mu2(mu2), xi(xi), levGran(a), pravGran(b)
	{
	}
	void calculate_test(int n)
	{
		x.clear();
		v.clear();
		v2.clear();
		diff.clear();

		vector<double> a;
		vector<double> d;
		vector<double> phi;

		double h = (pravGran - levGran) / n;
		for (int i = 0; i < n; i++)
		{
			a.push_back(1 / ((1 / k_test(levGran + h * i, xi) + k_test(levGran + h * (i + 1.0), xi)) / 2));
		}

		for (int i = 0; i < n - 1; i++)
		{
			phi.push_back((f_test(levGran + h * (i + 0.5), xi) + f_test(levGran + h * (i + 1.5), xi)) / 2);
			d.push_back((q_test(levGran + h * (i + 0.5), xi) + q_test(levGran + h * (i + 1.5), xi)) / 2);
		}

		vector<double> A;
		vector<double> C;
		vector<double> B;
		for (int i = 0; i < n - 1; i++)
		{
			A.push_back(a[i] / pow(h, 2));
			C.push_back((a[i] + a[i + 1]) / pow(h, 2) + d[i]);
			B.push_back(a[i + 1] / pow(h, 2));
		}

		for (int i = 0; i < n + 1; i++)
		{
			x.push_back(levGran + h * i);
		}
		v = solveMatrix(A, C, B, phi, mu1, mu2, n);
		v2 = calc_true_sol(levGran, pravGran, xi, n);
		for (int i = 0; i < n + 1; i++)
		{
			diff.push_back(std::fabs(v[i] - v2[i]));
		}
	}
	void calculate_main(int n)
	{
		x.clear();
		v.clear();
		v2.clear();
		diff.clear();

		vector<double> a1;
		vector<double> d1;
		vector<double> phi1;
		vector<double> a2;
		vector<double> d2;
		vector<double> phi2;

		double h1 = (pravGran - levGran) / n;
		double h2 = (pravGran - levGran) / (2.0*n);
		for (int i = 0; i < n; i++)
		{
			a1.push_back(1 / ((1 / k_main(levGran + h1 * i, xi) + k_main(levGran + h1 * (i + 1.0), xi)) / 2));
		}
		for (int i = 0; i < 2*n; i++)
		{
			a2.push_back(1 / ((1 / k_main(levGran + h2 * i, xi) + k_main(levGran + h2 * (i + 1.0), xi)) / 2));
		}

		for (int i = 0; i < n - 1; i++)
		{
			phi1.push_back((f_main(levGran + h1 * (i + 0.5), xi) + f_main(levGran + h1 * (i + 1.5), xi)) / 2);
			d1.push_back((q_main(levGran + h1 * (i + 0.5), xi) + q_main(levGran + h1 * (i + 1.5), xi)) / 2);
		}
		for (int i = 0; i < 2*n - 1; i++)
		{
			phi2.push_back((f_main(levGran + h2 * (i + 0.5), xi) + f_main(levGran + h2 * (i + 1.5), xi)) / 2);
			d2.push_back((q_main(levGran + h2 * (i + 0.5), xi) + q_main(levGran + h2 * (i + 1.5), xi)) / 2);
		}

		vector<double> A1;
		vector<double> C1;
		vector<double> B1;
		for (int i = 0; i < n - 1; i++)
		{
			A1.push_back(a1[i] / pow(h1, 2));
			C1.push_back((a1[i] + a1[i + 1]) / pow(h1, 2) + d1[i]);
			B1.push_back(a1[i + 1] / pow(h1, 2));
		}

		vector<double> A2;
		vector<double> C2;
		vector<double> B2;
		for (int i = 0; i < 2*n - 1; i++)
		{
			A2.push_back(a2[i] / pow(h2, 2));
			C2.push_back((a2[i] + a2[i + 1]) / pow(h2, 2) + d2[i]);
			B2.push_back(a2[i + 1] / pow(h2, 2));
		}

		for (int i = 0; i < n + 1; i++)
		{
			x.push_back(levGran + h1 * i);
		}
		v = solveMatrix(A1, C1, B1, phi1, mu1, mu2, n);
		vector<double> temp_v2 = solveMatrix(A2, C2, B2, phi2, mu1, mu2, 2*n);
		for (int i = 0; i < n+1; i++)
		{
			v2.push_back(temp_v2[2 * i]);
		}
		for (int i = 0; i < n + 1; i++)
		{
			diff.push_back(std::fabs(v[i] - v2[i]));
		}
	}
};