#pragma once
#include <vector>
#include <stdexcept>
#include <format>
#include <iostream>

#include "SparseMatrix.h"
#include "LU.h"

namespace Vec {
   inline double Scalar(const std::vector<double>& l, const std::vector<double>& r);

   // l or r may be similar vectors to ans
   inline void Mult(const std::vector<double>& l, const std::vector<double>& r, std::vector<double>& ans);
   inline std::vector<double> Mult(const std::vector<double>& l, const std::vector<double>& r);
}

namespace IterSolvers {
   extern double minEps;
   extern size_t maxIter;
   extern bool globalDebugOutput;

   extern std::vector<double>* _tmp1, * _tmp2,
      * _tmp3, * _tmp4, * _tmp5, * _tmp6;
   extern LU* _lu_mat;

   inline void VecInit(std::vector<double>*& vec, size_t size);

   namespace MSG_Assimetric {
      void Init_Default(size_t size);

      size_t Default(const SparseMatrix& A, const std::vector<double>& f, std::vector<double>& x, double& eps, bool debugOutput = globalDebugOutput);

      void Init_DiagPrecond(size_t size);

      size_t DiagPrecond(const SparseMatrix& A, const std::vector<double>& f, std::vector<double>& x, double& eps, bool debugOutput = globalDebugOutput);

      void Init_LuPrecond(size_t diSize, const SparseMatrix& A);

      size_t LuPrecond(const SparseMatrix& A, const std::vector<double>& f, std::vector<double>& x, double& eps, bool debugOutput = globalDebugOutput);
   }

   namespace LOS {
      extern size_t resetIter;

      void Init_Default(size_t size);

      size_t Default(const SparseMatrix& A, const std::vector<double>& f, std::vector<double>& x, double& eps, bool debugOutput = globalDebugOutput);

      void Init_DiagPrecond(size_t size);

      size_t DiagPrecond(const SparseMatrix& A, const std::vector<double>& f, std::vector<double>& x, double& eps, bool debugOutput = globalDebugOutput);

      void Init_LuPrecond(size_t diSize, const SparseMatrix& A);

      size_t LuPrecond(const SparseMatrix& A, const std::vector<double>& f, std::vector<double>& x, double& eps, bool debugOutput = globalDebugOutput);
   }

   void Destruct() noexcept;
};
