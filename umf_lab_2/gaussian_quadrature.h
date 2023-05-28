#pragma once
#include <functional>
#include <cmath>

namespace Gaussian_2p {
   constexpr int numPoint = 2;

   constexpr double weights[] = {1.0, 1.0};
   constexpr double points[] = {
      -0.577350269189625764509148780,    // -1 / sqrt(3)
      0.577350269189625764509148780,     //  1 / sqrt(3)
   };

   /// <summary>
   /// Class for compute integral of one-dimentional function, like \int_0^2 {2x + x dx}
   /// </summary>
   class OneDimentionSolver {
   public:
      double from = 0.0;
      double to = 0.0;

      std::function<double(double)> computeFunc;

   public:
      OneDimentionSolver(double from, double to, std::function<double(double)> func) :
         from(from), to(to), computeFunc(func) {}

      static OneDimentionSolver withStep(double from, double step, std::function<double(double)> func) {
         return OneDimentionSolver(from, from + step, func);
      }

   public:
      /// <summary>
      /// Function to compute quadrature with other range. Doesn't override previous ranges
      /// </summary>
      double computeWithRange(double from, double to) {
         const double coefs[] = {
            (to - from) / 2.0,
            (to + from) / 2.0,
         };

         double result = 0.0;
         for (auto i = 0; i < numPoint; i++)
         {
            result += weights[i] * computeFunc(coefs[0] * points[i] + coefs[1]);
         }
         return coefs[0] * result;
      }

      /// <summary>
      /// Function to compute quadrature by this equation:
      /// (to-from) / 2 * (w1 * f((to-from)/2 * x1 + (to+from)/2 ) + w2 * f((to-from)/2 * x1 + (to+from)/2 )
      /// </summary>
      /// <returns>result of computation (integral from [from] to [to])</returns>
      double compute() {
         return computeWithRange(from, to);
      }

      /// <summary>
      /// Function to compute quadrature with other range with step. Doesn't override previous ranges
      /// </summary>
      double computeWithStep(double from, double step) {
         return computeWithRange(from, from + step);
      }
   };

   /// <summary>
   /// Class for compute two-dimentional functions like int_xfrom^xTo {int_yFrom^yTo {x^2 + y^2 dy} dx}
   /// </summary>
   class TwoDimentionalSolver {
   public:
      double xFrom = 0.0;
      double xTo = 0.0;
      double yFrom = 0.0;
      double yTo = 0.0;

      std::function<double(double, double)> computeFunc;

   public:
      TwoDimentionalSolver(double xFrom, double xTo, double yFrom, double yTo, std::function<double(double, double)> computeFunc) :
         xFrom(xFrom), xTo(xTo), yFrom(yFrom), yTo(yTo), computeFunc(computeFunc) {}

      static TwoDimentionalSolver withStep(double xFrom, double xStep, double yFrom, double yStep, std::function<double(double, double)> computeFunc) {
         return TwoDimentionalSolver(xFrom, xFrom + xStep, yFrom, yFrom + yStep, computeFunc);
      }

   public:
      /// <summary>
      /// Function that return computed integral with fixed range, doesn't change in-object range
      /// </summary>
      /// <returns>Result of computation</returns>
      double computeWithRange(double xFrom, double xTo, double yFrom, double yTo) {
         const double x_coefs[] = {
            (xTo - xFrom) / 2.0,
            (xTo + xFrom) / 2.0,
         };
         const double y_coefs[] = {
            (yTo - yFrom) / 2.0,
            (yTo + yFrom) / 2.0,
         };

         double result = 0.0;
         for (auto i = 0; i < numPoint; i++)
         {
            for (auto j = 0; j < numPoint; j++)
            {
               result += weights[i] * weights[j]
                  * computeFunc(x_coefs[0] * points[i] + x_coefs[1], y_coefs[0] * points[j] + y_coefs[1]);
            }
         }
         return x_coefs[0] * y_coefs[0] * result;
      }

      /// <summary>
      /// Function that return computed integral with fixed in-object range
      /// </summary>
      /// <returns>Result of computation</returns>
      double compute() {
         return computeWithRange(xFrom, xTo, yFrom, yTo);
      }

      /// <summary>
      /// Function that return computed integral with fixed range, doesn't change in-object range
      /// </summary>
      /// <returns>Result of computation</returns>
      double computeWithStep(double xFrom, double xStep, double yFrom, double yStep) {
         return computeWithRange(xFrom, xFrom + xStep, yFrom, yFrom + yStep);
      }
   };
}

namespace Gaussian_4p {
   constexpr int numPoint = 4;

   constexpr double weights[] = {
      0.347854845137453857373063949222,
      0.652145154862546142626936050778,
      0.652145154862546142626936050778,
      0.347854845137453857373063949222
   };
   constexpr double points[] = {
      -0.86113631159405257522394648889281,
      -0.33998104358485626480266575910324,
      0.33998104358485626480266575910324,
      0.86113631159405257522394648889281,
   };

   /// <summary>
   /// Class for compute integral of one-dimentional function, like \int_0^2 {2x + x dx}
   /// </summary>
   class OneDimentionSolver {
   public:
      double from = 0.0;
      double to = 0.0;

      std::function<double(double)> computeFunc;

   public:
      OneDimentionSolver(double from, double to, std::function<double(double)> func) :
         from(from), to(to), computeFunc(func) {}

      static OneDimentionSolver withStep(double from, double step, std::function<double(double)> func) {
         return OneDimentionSolver(from, from + step, func);
      }

   public:
      /// <summary>
      /// Function to compute quadrature with other range. Doesn't override previous ranges
      /// </summary>
      double computeWithRange(double from, double to) {
         const double coefs[] = {
            (to - from) / 2.0,
            (to + from) / 2.0,
         };

         double result = 0.0;
         for (auto i = 0; i < numPoint; i++)
         {
            result += weights[i] * computeFunc(coefs[0] * points[i] + coefs[1]);
         }
         return coefs[0] * result;
      }

      /// <summary>
      /// Function to compute quadrature by this equation:
      /// (to-from) / 2 * (w1 * f((to-from)/2 * x1 + (to+from)/2 ) + w2 * f((to-from)/2 * x1 + (to+from)/2 )
      /// </summary>
      /// <returns>result of computation (integral from [from] to [to])</returns>
      double compute() {
         return computeWithRange(from, to);
      }

      /// <summary>
      /// Function to compute quadrature with other range with step. Doesn't override previous ranges
      /// </summary>
      double computeWithStep(double from, double step) {
         return computeWithRange(from, from + step);
      }
   };

   /// <summary>
   /// Class for compute two-dimentional functions like int_xfrom^xTo {int_yFrom^yTo {x^2 + y^2 dy} dx}
   /// </summary>
   class TwoDimentionalSolver {
   public:
      double xFrom = 0.0;
      double xTo = 0.0;
      double yFrom = 0.0;
      double yTo = 0.0;

      std::function<double(double, double)> computeFunc;

   public:
      TwoDimentionalSolver(double xFrom, double xTo, double yFrom, double yTo, std::function<double(double, double)> computeFunc) :
         xFrom(xFrom), xTo(xTo), yFrom(yFrom), yTo(yTo), computeFunc(computeFunc) {}

      static TwoDimentionalSolver withStep(double xFrom, double xStep, double yFrom, double yStep, std::function<double(double, double)> computeFunc) {
         return TwoDimentionalSolver(xFrom, xFrom + xStep, yFrom, yFrom + yStep, computeFunc);
      }

   public:
      /// <summary>
      /// Function that return computed integral with fixed range, doesn't change in-object range
      /// </summary>
      /// <returns>Result of computation</returns>
      double computeWithRange(double xFrom, double xTo, double yFrom, double yTo) {
         const double x_coefs[] = {
            (xTo - xFrom) / 2.0,
            (xTo + xFrom) / 2.0,
         };
         const double y_coefs[] = {
            (yTo - yFrom) / 2.0,
            (yTo + yFrom) / 2.0,
         };

         double result = 0.0;
         for (auto i = 0; i < numPoint; i++)
         {
            for (auto j = 0; j < numPoint; j++)
            {
               result += weights[i] * weights[j]
                  * computeFunc(x_coefs[0] * points[i] + x_coefs[1], y_coefs[0] * points[j] + y_coefs[1]);
            }
         }
         return x_coefs[0] * y_coefs[0] * result;
      }

      /// <summary>
      /// Function that return computed integral with fixed in-object range
      /// </summary>
      /// <returns>Result of computation</returns>
      double compute() {
         return computeWithRange(xFrom, xTo, yFrom, yTo);
      }

      /// <summary>
      /// Function that return computed integral with fixed range, doesn't change in-object range
      /// </summary>
      /// <returns>Result of computation</returns>
      double computeWithStep(double xFrom, double xStep, double yFrom, double yStep) {
         return computeWithRange(xFrom, xFrom + xStep, yFrom, yFrom + yStep);
      }
   };
}