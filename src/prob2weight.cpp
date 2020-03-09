#include <stdio.h>
#include <vector>

using namespace std;

#include "prob2weight.h"

void test_prob2weight()
{
	const int K = 141235;
	double a = 4, ta = 4;

	Prob2Weight<> p2w(K, pow(a, ta));

	printf("K = %d\n", p2w.K);
	printf("W = %d\n", p2w.W);
	printf("WE = %d\n", p2w.WE);
	printf("A = %f\n", p2w.A);
	printf("p0 = %e\n", p2w.p0);
	printf("L = %f\n", p2w.L);
	printf("G0 = %f\n", p2w.G0);
	printf("G1 = %f\n", p2w.G1);
	printf("scatter type = %s\n", Smooth == p2w.ScatterType ? "smooth" : "accurate boundary");
	printf("D = %f\n", p2w.D);
	printf("tk = %f\n", p2w.tk);

	printf("\n");
	printf("prob %d -> origin weight %f\n", 0, p2w.g(0));
	printf("prob %d -> origin weight %f\n", 1, p2w.g(1));
	printf("prob %d -> weight %f\n", 0, p2w.f(0));
	printf("prob %d -> weight %f\n", 1, p2w.f(1));

	printf("\n");
	vector<double> probs;
	for (int i = -1*ta; i <= ta; i++)
		probs.push_back(pow(a, i) * p2w.p0);
	for (auto p : probs)
		printf("prob %e -> weight %f\n", p, p2w.f(p));
	printf("\n");
	vector<double> prob_multis;
	for (int i = 1; i <= 10; i++)
		prob_multis.push_back(i*0.1);
	for (auto m : prob_multis)
		printf("prob multi %f -> weight %f\n", m, p2w.f_m(m));
}

int main()
{
	printf("-------- test --------\n");
	test_prob2weight();

	printf("-------- demo --------\n");
	int K = 141235;
	Prob2Weight<> p2w(K);
	vector<double> probs = {0, 0.01/K, 0.1/K, 1.0/K, 10.0/K, 100.0/K, 1.0/2, 1};
	for (auto p : probs)
		printf("prob %e -> weight %d\n", p, p2w.Weight(p));
	vector<double> prob_multis = {0.25, 0.5, 0.75, 1};
	printf("\n");
	for (auto m : prob_multis)
		printf("prob multi %f -> weight %d\n", m, p2w.MultipleWeight(m));

	return 0;
}