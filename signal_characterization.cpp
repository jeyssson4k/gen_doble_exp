/*
* Double exponential lightning-like signal characterization
* @Author Jeysson4K
*/

#include <stdio.h>
#include <math.h>
#include <stdexcept>
#include <iostream>
#include <vector>

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
#define CAPACITOR_A 1e-6
#define CAPACITOR_B 3.3e-6

#define EX_T1 1e-3
#define EX_T2 10e-3
#define EX_ALPHA 0.05
#define EX_BETA 1.0

void tau2_equation(double &tau2, const double T2, const double alpha, const double A, const double B);
double f(const double tau, const double T);
std::vector<double> cascade_filter_params(const double tau1, const double tau2, const double B, const double cb, const double ca);

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
    T1 = T1_MS;
    T2 = T2_MS;
    alpha = 0.5;

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


    printf("\nParametros circuito pasabanda en cascada\n");
    std::vector<double> pr = cascade_filter_params(tau1, tau2, B, CAPACITOR_B, CAPACITOR_A);
    printf("{ Cb: %.2e,\tRc: %.2f,\tRd: %.2f,\tCa: %.2e,\tRb: %.2f,\tRa: %.2f }\n", pr[0], pr[1], pr[2], pr[3], pr[4], pr[5]);
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

std::vector<double> cascade_filter_params(const double tau1, const double tau2, const double B, const double cb, const double ca)
{
    std::vector<double> params{cb, .0, .0, ca, .0, .0};

    params[1] = tau2/cb;
    params[2] = B*params[1];
    params[4] = tau1/ca;
    params[5] = (params[2]*(tau1*tau2))/((params[1]*ca*B)*(tau2-tau1));

    return params; 
}
