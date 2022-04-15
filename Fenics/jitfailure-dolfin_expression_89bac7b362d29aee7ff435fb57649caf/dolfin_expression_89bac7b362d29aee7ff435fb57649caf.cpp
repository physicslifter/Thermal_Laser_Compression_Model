
// Based on https://gcc.gnu.org/wiki/Visibility
#if defined _WIN32 || defined __CYGWIN__
    #ifdef __GNUC__
        #define DLL_EXPORT __attribute__ ((dllexport))
    #else
        #define DLL_EXPORT __declspec(dllexport)
    #endif
#else
    #define DLL_EXPORT __attribute__ ((visibility ("default")))
#endif

#include <dolfin/function/Expression.h>
#include <dolfin/math/basic.h>
#include <Eigen/Dense>


// cmath functions
using std::cos;
using std::sin;
using std::tan;
using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::cosh;
using std::sinh;
using std::tanh;
using std::exp;
using std::frexp;
using std::ldexp;
using std::log;
using std::log10;
using std::modf;
using std::pow;
using std::sqrt;
using std::ceil;
using std::fabs;
using std::floor;
using std::fmod;
using std::max;
using std::min;

const double pi = DOLFIN_PI;


namespace dolfin
{
  class dolfin_expression_89bac7b362d29aee7ff435fb57649caf : public Expression
  {
     public:
       double tol;
double k_1;
std::shared_ptr<dolfin::GenericFunction> generic_function_u;
double a;
double b;
double c;
double l;


       dolfin_expression_89bac7b362d29aee7ff435fb57649caf()
       {
            
       }

       void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
       {
          double u;
            generic_function_u->eval(Eigen::Map<Eigen::Matrix<double, 1, 1>>(&u), x);
          values[0] = x[0] <= l + tol ? b+a*u+c*d : k_1;

       }

       void set_property(std::string name, double _value) override
       {
          if (name == "tol") { tol = _value; return; }          if (name == "k_1") { k_1 = _value; return; }          if (name == "a") { a = _value; return; }          if (name == "b") { b = _value; return; }          if (name == "c") { c = _value; return; }          if (name == "l") { l = _value; return; }
       throw std::runtime_error("No such property");
       }

       double get_property(std::string name) const override
       {
          if (name == "tol") return tol;          if (name == "k_1") return k_1;          if (name == "a") return a;          if (name == "b") return b;          if (name == "c") return c;          if (name == "l") return l;
       throw std::runtime_error("No such property");
       return 0.0;
       }

       void set_generic_function(std::string name, std::shared_ptr<dolfin::GenericFunction> _value) override
       {
          if (name == "u") { generic_function_u = _value; return; }
       throw std::runtime_error("No such property");
       }

       std::shared_ptr<dolfin::GenericFunction> get_generic_function(std::string name) const override
       {
          if (name == "u") return generic_function_u;
       throw std::runtime_error("No such property");
       }

  };
}

extern "C" DLL_EXPORT dolfin::Expression * create_dolfin_expression_89bac7b362d29aee7ff435fb57649caf()
{
  return new dolfin::dolfin_expression_89bac7b362d29aee7ff435fb57649caf;
}

