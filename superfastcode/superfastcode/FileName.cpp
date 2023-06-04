#include <Windows.h>
#include <cmath>
#include <Python.h>
#include <iostream>
#include <cmath>
# define M_PI           3.14159265358979323846  /* pi */
double K = 95;
double T = 0.25;
double S = 100;
double r = 0.075;
double market_price = 10;
bool type = 0;
double norm_pdf(const double x) {
    return (1.0 / (pow(2 * M_PI, 0.5))) * exp(-0.5 * x * x);
}

// An approximation to the cumulative distribution function
// for the standard normal distribution
// Note: This is a recursive function
double norm_cdf(const double x) {
    double k = 1.0 / (1.0 + 0.2316419 * x);
    double k_sum = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937 + k * (-1.821255978 + 1.330274429 * k))));

    if (x >= 0.0) {
        return (1.0 - (1.0 / (pow(2 * M_PI, 0.5))) * exp(-0.5 * x * x) * k_sum);
    }
    else {
        return 1.0 - norm_cdf(-x);
    }
}
double d_j(const int j, const double S, const double K, const double r, const double sigma, const double T)
{
    return (log(S / K) + (r + (pow(-1, j)) * 0.5 * sigma * sigma) * T) / (sigma * (pow(T, 0.5)));
}

double call_price(const double S, const double K, const double r, const double sigma, const double T)
{
    return (S * norm_cdf(d_j(2, S, K, r, sigma, T)) - K * exp(-r * T) * norm_cdf(d_j(1, S, K, r, sigma, T)));
}

double put_price(const double S, const double K, const double r, const double sigma, const double T)
{
    return (K * exp(-r * T) * norm_cdf((-1) * d_j(1, S, K, r, sigma, T)) - S * norm_cdf((-1) * d_j(2, S, K, r, sigma, T)));
}

double call_vega(const double S, const double K, const double r, const double sigma, const double T) {
    return S * sqrt(T) * norm_pdf(d_j(1, S, K, r, sigma, T));
}




double newton_raphson()


{
    double (*model)(const double S, const double K, const double r, const double sigma, const double T);

    if (type == 0) {
        model = &call_price;
    }
    else
    {
        model = &put_price;
    }

    double epsilon = 1e-7;
    double IV = 0.8;

    double model_price = model(S, K, r, IV, T);
    double vega;
    double diff;

    while (fabs(model_price - market_price) > epsilon) {

        model_price = model(S, K, r, IV, T);
        vega = call_vega(S, K, r, IV, T);
        diff = market_price - model_price;
        IV = IV + diff / vega;
        std::cout << IV << std::endl;

    }

    return IV;

}




PyObject* tanh_impl(PyObject* self, PyObject *o)
{
    PyArg_ParseTuple(o, "pddddd", &type, &market_price,&S, &K, &r, &T);


    /*
    

( Type S, K, r, IV, T);*/

    double result= newton_raphson();

    return PyFloat_FromDouble(result);
}





static PyMethodDef superfastcode_methods[] = {
    // The first property is the name exposed to Python, fast_tanh
    // The second is the C++ function with the implementation
    // METH_O means it takes a single PyObject argument
    { "fast_tanh", (PyCFunction)tanh_impl, METH_VARARGS, nullptr },

    // Terminate the array with an object containing nulls.
    { nullptr, nullptr, 0, nullptr }
};

static PyModuleDef superfastcode_module = {
    PyModuleDef_HEAD_INIT,
    "superfastcode",                        // Module name to use with Python import statements
    "Provides some functions, but faster",  // Module description
    0,
    superfastcode_methods                   // Structure that defines the methods of the module
};

PyMODINIT_FUNC PyInit_superfastcode() {
    return PyModule_Create(&superfastcode_module);
}
