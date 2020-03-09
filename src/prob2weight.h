#include <cmath>

using namespace std;

#define _A 256

enum ScatterType {
	Smooth = 0,
	AccurateBoundary,
};

template<int _W = 256, int _WE = 20, int _ScatterType = Smooth>
class Prob2Weight
{
public:
	Prob2Weight(int K, double A = _A, double p0_K = 1.0)
	: K(K), A(A), p0(p0_K/K)
	{
		L = (W - 2*WE) / (2*log(A));

		const double delta = 0.5;
		D = -1 * log(exp(delta/L) - 1);
		if ( Smooth == ScatterType )
			tk = ((W-2)/(2*L) - log(p0)) / log(K); // smooth
		else
			tk = (log(A) - log(p0) + D) / log(K); // accurate boundary

		_K_tk = pow(K, tk);

		_h_p0 = h(p0);

		_p0_A_ = p0/A;
		_p0_A = p0*A;
		G0 = g(_p0_A_);
		G1 = g(_p0_A);

		_W_1_G0_g_0_G0 = (W-1-G0)/(g(0)-G0);
		_1_G1_g_1_G1 = (1-G1)/(g(1)-G1);
	}

	int Weight(double p) const
	{ return int(f(p)+0.5); }

	// Weight(x*m) = Weight(x) + MultipleWeight(m)
	int MultipleWeight(double m) const
	{ return int(f_m(m)+0.5); }

private:
	static const int W = _W;
	static const int WE = _WE;
	static const int ScatterType = _ScatterType;
	
	inline double h(double x) const
	{ return log(x*_K_tk + 1); }
	inline double g(double x) const
	{ return (_h_p0 - h(x))*L + 0.5*W; }
	inline double f(double x) const
	{
		double g_x = g(x);
		if (x < _p0_A_)
			g_x = G0 + (g_x-G0)*_W_1_G0_g_0_G0;
		else
		if (_p0_A < x)
			g_x = G1 + (g_x-G1)*_1_G1_g_1_G1;
		return g_x;
	}
	inline double f_m(double x) const
	{ return -1*log(x)*L; }

	friend void test_prob2weight();

private:
	int K;
	double A, p0;
	double tk, D;
	double L;
	double G0, G1;

	double _K_tk, _h_p0;
	double _p0_A_, _p0_A;
	double _W_1_G0_g_0_G0, _1_G1_g_1_G1;
};
