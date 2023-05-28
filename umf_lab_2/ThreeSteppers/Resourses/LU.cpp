#include "../Headers/LU.h"

/// <summary>
/// Конструктор с резервированием памяти под разложение
/// </summary>
/// <param name="diSize"> - размер диагонали,</param>
/// <param name="luSize"> - размер массивов нижнего и верхнего треугольника</param>
LU::LU(size_t diSize, size_t luSize) {
   Resize(diSize, luSize);
}

/// <summary>
/// Конструктор с построением неполного LU(sq)-разложения по матрице mat
/// </summary>
/// <param name="mat"> - матрица, по которой построится LU-разложение, с привязкой этой матрицы к объекту</param>
LU::LU(const SparseMatrix& mat)
{
   MakeLuFor(mat);
}

/// <summary>
/// Разложить матрицу mat в неполное LU(sq) - разложение
/// </summary>
/// <param name="mat"> - матрица, которую требуется разложить. Она же будет использоваться для просмотра портрета матриц</param>
void LU::MakeLuFor(const SparseMatrix& mat) {
   parent = &mat;
   if (di.size() != mat.di.size())
      di.resize(mat.di.size());
   if (ggl.size() != mat.ggl.size())
      ggl.resize(mat.ggl.size());
   if (ggu.size() != mat.ggu.size())
      ggu.resize(mat.ggu.size());

   const auto& ig = mat.ig;
   const auto& jg = mat.jg;

   for (size_t i = 0; i < mat.Size(); i++)
   {
      double di_accum = 0;
      for (size_t j = ig[i]; j < ig[i + 1]; j++)
      {
         size_t k = ig[i];
         size_t v = ig[jg[j]];
         double ggl_accum = 0, ggu_accum = 0;
         while (k < j && v < ig[jg[j] + 1ll])
         {
            if (jg[k] > jg[v]) v++;
            else if (jg[k] < jg[v]) k++;
            else
            {
               ggl_accum += ggl[k] * ggu[v];
               ggu_accum += ggl[v] * ggu[k];
               k++;
               v++;
            }
         }
         ggl[j] = (mat.ggl[j] - ggl_accum) / di[jg[j]];
         ggu[j] = (mat.ggu[j] - ggu_accum) / di[jg[j]];

         di_accum += ggl[j] * ggu[j];
      }

      di[i] = sqrt(mat.di[i] - di_accum);
   }
}

void LU::Resize(size_t diSize, size_t luSize) {
   di.resize(diSize);
   ggl.resize(luSize);
   ggu.resize(luSize);
}

/// <summary>
/// Умножение нижней матрицы L на вектор vec. Выделяет память под вектор ответа, не меняет матрицу LU
/// </summary>
/// <param name="vec"> - вектор, на который будет происходить умножение матрицы L</param>
/// <returns>Вектор с результатом перемножения (выделяется в памяти)</returns>
std::vector<double> LU::LMultToVec(const std::vector<double>& vec) const
{
   std::vector<double> ans(vec.size());
   return LMultToVec(vec, ans);
}

/// <summary>
/// Умножение нижней матрицы L на вектор vec. Ответ записывается в вектор ans, не меняет матрицу LU
/// </summary>
/// <param name="vec"> - вектор, на который будет происходить умножение матрицы L;</param>
/// <param name="ans"> - вектор, куда запишется ответ без выделения памяти (должен отличаться от vec!)</param>
/// <returns>Ссылка на вектор ans</returns>
std::vector<double>& LU::LMultToVec(const std::vector<double>& vec, std::vector<double>& ans) const
{
   if (vec.size() != di.size()) throw std::runtime_error("Размеры матрицы и вектора не совпадают.");
   if (vec.size() != ans.size()) throw std::runtime_error("Размеры матрицы и результирующего вектора не совпадают.");

   for (uint16_t i = 0; i < ans.size(); i++)
   {
      // Умножаем диагональ
      ans[i] = di[i] * vec[i];

      // Умножаем нижний треугольник
      for (uint32_t j = parent->ig[i]; j < parent->ig[i + 1ll]; j++)
      {
         ans[i] += ggl[j] * vec[parent->jg[j]];
      }
   }

   return ans;
}

std::vector<double> LU::LTranspMultToVec(const std::vector<double>& vec) const
{
   std::vector<double> ans(vec.size());
   return LTranspMultToVec(vec, ans);
}

std::vector<double>& LU::LTranspMultToVec(const std::vector<double>& vec, std::vector<double>& ans) const
{
   if (vec.size() != di.size()) throw std::runtime_error("Размеры матрицы и вектора не совпадают.");
   if (vec.size() != ans.size()) throw std::runtime_error("Размеры матрицы и результирующего вектора не совпадают.");

   for (uint16_t i = 0; i < ans.size(); i++)
   {
      // Умножаем диагональ
      ans[i] = di[i] * vec[i];

      // Умножаем на верхний треугольник с данными нижнего
      for (uint32_t j = parent->ig[i]; j < parent->ig[i + 1ll]; j++)
      {
         ans[parent->jg[j]] += ggl[j] * vec[i];
      }
   }

   return ans;
}

/// <summary>
/// Умножение верхней матрицы U на вектор vec. Выделяет память под вектор ответа, не меняет матрицу LU
/// </summary>
/// <param name="vec"> - вектор, на который будет происходить умножение матрицы U</param>
/// <returns>Вектор с результатом перемножения (выделяется в памяти)</returns>
std::vector<double> LU::UMultToVec(const std::vector<double>& vec) const
{
   std::vector<double> ans(vec.size());
   return UMultToVec(vec, ans);
}

/// <summary>
/// Умножение верхней матрицы U на вектор vec. Ответ записывается в вектор ans, не меняет матрицу LU
/// </summary>
/// <param name="vec"> - вектор, на который будет происходить умножение матрицы U;</param>
/// <param name="ans"> - вектор, куда запишется ответ без выделения памяти (должен отличаться от vec!)</param>
/// <returns>Ссылка на вектор ans</returns>
std::vector<double>& LU::UMultToVec(const std::vector<double>& vec, std::vector<double>& ans) const
{
   if (vec.size() != di.size()) throw std::runtime_error("Размеры матрицы и вектора не совпадают.");
   if (vec.size() != ans.size()) throw std::runtime_error("Размеры матрицы и результирующего вектора не совпадают.");

   for (uint16_t i = 0; i < ans.size(); i++)
   {
      // Умножаем диагональ
      ans[i] = di[i] * vec[i];

      // Умножаем верхний треугольник
      for (uint32_t j = parent->ig[i]; j < parent->ig[i + 1ll]; j++)
      {
         ans[parent->jg[j]] += ggu[j] * vec[i];
      }
   }

   return ans;
}

std::vector<double> LU::UTranspMultToVec(const std::vector<double>& vec) const
{
   std::vector<double> ans(vec.size());
   return UTranspMultToVec(vec, ans);
}

std::vector<double>& LU::UTranspMultToVec(const std::vector<double>& vec, std::vector<double>& ans) const
{
   if (vec.size() != di.size()) throw std::runtime_error("Размеры матрицы и вектора не совпадают.");
   if (vec.size() != ans.size()) throw std::runtime_error("Размеры матрицы и результирующего вектора не совпадают.");

   for (uint16_t i = 0; i < ans.size(); i++)
   {
      // Умножаем диагональ
      ans[i] = di[i] * vec[i];

      // Умножаем нижний треугольник с данными верхнего треугольника
      for (uint32_t j = parent->ig[i]; j < parent->ig[i + 1ll]; j++)
      {
         ans[i] += ggu[j] * vec[parent->jg[j]];
      }
   }

   return ans;
}

/// <summary>
/// Решение слау вида Lx = right. Не выделяет память под вектор x, не меняет матрицы LU
/// </summary>
/// <param name="right"> - вектор правой части уравнения;</param>
/// <param name="x"> - вектор, куда будет записан ответ. Должен быть с уже выделенной памятью. Должен отличаться от right!</param>
/// <returns>ссылка на вектор x</returns>
std::vector<double>& LU::LSlauSolve(const std::vector<double>& right, std::vector<double>& x) const
{
   if (right.size() != di.size()) throw std::runtime_error("Размеры матрицы и вектора не совпадают.");
   if (right.size() != x.size()) throw std::runtime_error("Размеры матрицы и результирующего вектора не совпадают.");

   size_t size = x.size();
   for (size_t i = 0; i < size; i++)
   {
      x[i] = 0;
      for (size_t j = parent->ig[i]; j < parent->ig[i + 1]; j++)
      {
         x[i] += x[parent->jg[j]] * ggl[j];
      }
      x[i] = (right[i] - x[i]) / di[i];
   }

   return x;
}

/// <summary>
/// Решение слау вида Lx = right. Выделяет память под вектор x, не меняет матрицы LU
/// </summary>
/// <param name="right"> - вектор правой части уравнения;</param>
/// <returns>Полученный вектор x</returns>
std::vector<double> LU::LSlauSolve(const std::vector<double>& right) const
{
   std::vector<double> x(right.size());
   return LSlauSolve(right, x);
}

std::vector<double>& LU::LTranspSlauSolve(const std::vector<double>& right, std::vector<double>& x) const
{
   if (right.size() != di.size()) throw std::runtime_error("Размеры матрицы и вектора не совпадают.");
   if (right.size() != x.size()) throw std::runtime_error("Размеры матрицы и результирующего вектора не совпадают.");

   size_t size = x.size();
   for (size_t i = 0; i < size; i++)
      x[i] = 0;

   for (size_t it = 0, i = size - 1; it < size; it++, i--)
   {
      x[i] = (right[i] - x[i]) / di[i];
      for (size_t j = parent->ig[i]; j < parent->ig[i + 1]; j++)
      {
         x[parent->jg[j]] += ggl[j] * x[i];
      }
   }

   return x;
}

std::vector<double> LU::LTranspSlauSolve(const std::vector<double>& right) const
{
   std::vector<double> x(right.size());
   return LTranspSlauSolve(right, x);
}

/// <summary>
/// Решение слау вида Ux = right. Не выделяет память под вектор x, не меняет матрицы LU
/// </summary>
/// <param name="right"> - вектор правой части уравнения;</param>
/// <param name="x"> - вектор, куда будет записан ответ. Должен быть с уже выделенной памятью. Должен отличаться от right!</param>
/// <returns>ссылка на вектор x</returns>
std::vector<double>& LU::USlauSolve(const std::vector<double>& right, std::vector<double>& x) const
{
   if (right.size() != di.size()) throw std::runtime_error("Размеры матрицы и вектора не совпадают.");
   if (right.size() != x.size()) throw std::runtime_error("Размеры матрицы и результирующего вектора не совпадают.");

   size_t size = x.size();
   for (size_t i = 0; i < size; i++)
      x[i] = 0;

   for (size_t it = 0, i = size - 1; it < size; it++, i--)
   {
      x[i] = (right[i] - x[i]) / di[i];
      for (size_t j = parent->ig[i]; j < parent->ig[i + 1]; j++)
      {
         x[parent->jg[j]] += ggu[j] * x[i];
      }
   }

   return x;
}

/// <summary>
/// Решение слау вида Ux = right. Выделяет память под вектор x, не меняет матрицы LU
/// </summary>
/// <param name="right"> - вектор правой части уравнения;</param>
/// <returns>Полученный вектор x</returns>
std::vector<double> LU::USlauSolve(const std::vector<double>& right) const
{
   std::vector<double> x(right.size());
   return USlauSolve(right, x);
}

std::vector<double>& LU::UTranspSlauSolve(const std::vector<double>& right, std::vector<double>& x) const
{
   if (right.size() != di.size()) throw std::runtime_error("Размеры матрицы и вектора не совпадают.");
   if (right.size() != x.size()) throw std::runtime_error("Размеры матрицы и результирующего вектора не совпадают.");

   size_t size = x.size();
   for (size_t i = 0; i < size; i++)
   {
      x[i] = 0;
      for (size_t j = parent->ig[i]; j < parent->ig[i + 1]; j++)
      {
         x[i] += x[parent->jg[j]] * ggu[j];
      }
      x[i] = (right[i] - x[i]) / di[i];
   }

   return x;
}

std::vector<double> LU::UTranspSlauSolve(const std::vector<double>& right) const
{
   std::vector<double> x(right.size());
   return UTranspSlauSolve(right, x);
}