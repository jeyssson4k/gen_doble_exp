/*
* Double exponential lightning-like signal characterization
* @Author Jeysson4K
*/

#include <stdio.h>
#include <math.h>
#include <stdexcept>
#include <iostream>

#define PEAK 1.0
#define T1_MS 0.5e-3
#define T1_US 0.5e-6
#define T2_MS 5e-3
#define T2_US 5e-6
#define DEFAULT_ALPHA 0.1
#define MAX_APPROX_IT 1000
#define MAX_IT 100
#define TOL 1e-300
#define DEFAULT 0.0
#define TEMP 1e+6
#define S 0.0

#define EX_T1 1e-3
#define EX_T2 10e-3
#define EX_ALPHA 0.05
#define EX_BETA 1.0

void tau2_equation(double &tau2, const double T2, const double alpha, const double A, const double B);
double f(const double tau, const double T);

int main(void)
{

    // params
    double
        beta,  // relationship A/B
        sbeta, // relationship to control approximation
        A,     // desired signal peak
        B,     // constant B[x1(t) - x2(t)]
        tau1,  // time of x1(t)
        tau2,  // time of x2(t)
        T1,    // period of x1(t)
        T2,    // period of x2(t)
        alpha; // control variable between (0,1)

    // initial conditions
    beta = EX_BETA;
    sbeta = DEFAULT;
    A = PEAK;
    B = DEFAULT;
    tau1 = DEFAULT;
    tau2 = DEFAULT;
    T1 = T1_US;
    T2 = T2_US;
    alpha = DEFAULT_ALPHA;

    // program-control varaibles
    bool t = false;
    int k = 1;

    // start computing tau1, tau2 and B
    printf("it\ttau1\t tau2\t B\t beta\n");
    do
    {
        if (t)
        {
            beta = sbeta;
        }

        // computing B
        B = A / beta;

        // computing tau2
        tau2_equation(tau2, T2, alpha, A, B);

        // computing tau1
        double tau_err = TEMP;
        double F_2 = f(tau2, T1);
        int N = MAX_APPROX_IT;
        for (int i = 2; i < N; ++i)
        {
            double target = tau2 / (S + i);
            double F_i = f(target, T1);

            double fx = (F_i - F_2) * (F_i - F_2);

            if (fx <= tau_err)
            {
                tau1 = target;
                tau_err = fx;
            }
        }

        double x = S - (T1 / tau2);
        double y = S - (T1 / tau1);
        sbeta = std::exp(x) - std::exp(y);
        printf("%i\t%.4e\t%.4e\t%.4e\t%.4e\n", k, tau1, tau2, B, beta);
        t = true;
        k++;
    } while (std::pow((beta - sbeta), 2) > TOL and k < MAX_IT);

    return EXIT_SUCCESS;
}

void tau2_equation(double &tau2, const double T2, const double alpha, const double A, const double B)
{
    double x = (alpha * A) / B;
    double y = log(x);
    tau2 = S - (T2 / y);
    return;
}

double f(const double tau, const double T)
{
    double x = S - (T / tau);
    return std::exp(x) / tau;
}
