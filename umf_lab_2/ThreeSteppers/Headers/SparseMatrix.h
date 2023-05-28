#pragma once
#include <cmath>
#include <format>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

std::vector<double> ReadVecFromFile(size_t size, const std::string& path);

/// <summary>
/// Класс объектов матриц, хранящихся в разреженном строчно-столбцовом виде
/// <para> Точность хранения элементов - double </para>
/// </summary>
class SparseMatrix {
   // Переменные матрицы
  public:
   /// <summary>
   /// Массив индексов строк/столбцов, вида 0, 0, 0 + k2, ..., 0+k2+...+kn, где ki - число элементов в i cтроке/столбце
   /// <para> Помимо этого первый элемент i строки можно найти как ggl[ig[i]] </para>
   /// <para> Пример массива для матрицы 3х3: </para>
   ///
   /// <para> Матрица: </para>
   /// <para> | 1 2 0 | </para>
   /// <para> | 3 8 1 | </para>
   /// <para> | 0 2 4 | </para>
   /// <para> ig: { 0, 0, 1, 2 } </para>
   /// </summary>
   std::vector<uint32_t> ig;

   /// <summary>
   /// Массив индексов столбцов/строк элементов (ставит индекс в соответствие элементу)
   /// <para> Пример массива для матрицы 3х3: </para>
   ///
   /// <para> Матрица: </para>
   /// <para> | 1 2 0 | </para>
   /// <para> | 3 8 1 | </para>
   /// <para> | 0 2 4 | </para>
   /// <para> jg: { 0, 1 } </para>
   /// </summary>
   std::vector<uint16_t> jg;

   /// <summary>
   /// Массив элементов нижнего треугольника матрицы
   /// <para> Пример массива для матрицы 3х3: </para>
   ///
   /// <para> Матрица: </para>
   /// <para> | 1 2 0 | </para>
   /// <para> | 3 8 1 | </para>
   /// <para> | 0 2 4 | </para>
   /// <para> ggl: { 3, 2 } </para>
   /// </summary>
   std::vector<double> ggl;

   /// <summary>
   /// Массив элементов верхнего треугольника матрицы
   /// <para> Пример массива для матрицы 3х3: </para>
   ///
   /// <para> Матрица: </para>
   /// <para> | 1 2 0 | </para>
   /// <para> | 3 8 1 | </para>
   /// <para> | 0 2 4 | </para>
   /// <para> ggu: { 2, 1 } </para>
   /// </summary>
   std::vector<double> ggu;

   /// <summary>
   /// Массив элементов диагонали матрицы
   /// <para> Пример массива для матрицы 3х3: </para>
   ///
   /// <para> Матрица: </para>
   /// <para> | 1 2 0 | </para>
   /// <para> | 3 8 1 | </para>
   /// <para> | 0 2 4 | </para>
   /// <para> di: { 1, 8, 4 } </para>
   /// </summary>
   std::vector<double> di;

   // Методы матрицы
  public:
   size_t Size() const;

   //Сложение разреженных матриц одной формы
   SparseMatrix operator+ (const SparseMatrix& other) const;

   /// <summary>
   /// Умножение матрицы на вектор
   /// </summary>
   std::vector<double> MultToVec(const std::vector<double>& right) const;
   std::vector<double>& MultToVec(const std::vector<double>& right, std::vector<double>& result) const;
   std::vector<double> operator*(const std::vector<double>& right) const;

   // Умножение матрицы на число
   SparseMatrix multToScalar(double scalar) const;

   // Умножение матрицы на число
   SparseMatrix operator*(double scalar) const { return multToScalar(scalar); }

   // Умножение матрицы на число
   friend SparseMatrix operator*(double scalar, const SparseMatrix& right) { return right * scalar; }

   /// <summary>
   /// Умножение транспонированной матрицы на вектор
   /// </summary>
   std::vector<double> TranspMultToVec(const std::vector<double>& right) const;
   std::vector<double>& TranspMultToVec(const std::vector<double>& right, std::vector<double>& result) const;

   SparseMatrix& operator=(SparseMatrix&& right) noexcept;

   double val(uint16_t row, uint16_t column);

   std::string toStringAsDense() {
      std::string out = "[ ";
      auto size = Size();

      for (auto i = 0; i < size; i++) {
         if (i != 0) out += "  ";
         out += "[ ";
         for (auto j = 0; j < size; j++) {
            out += std::format("{: 15.5f}", val(i, j));  // std::to_string(val(i, j));
            if (j + 1ll < size) out += ", ";
         }
         out += " ]";
         if (i + 1ll < size) out += "\n";
      }
      out += " ]";

      return out;
   }

   // Конструкторы матрицы
  public:
   SparseMatrix();

   // Конструктор копирования
   SparseMatrix(const SparseMatrix& right);

   // Конструктор перемещения (нужен для метода ReadFromFiles)
   SparseMatrix(SparseMatrix&& right) noexcept;

   // Статические методы матрицы
  public:
   // Конструктор копирования формы матрицы
   static SparseMatrix copyShape(const SparseMatrix& other);

   static SparseMatrix ReadFromFiles(uint16_t matrixSize, const std::string& igP, const std::string& jgP,
                                     const std::string& gglP, const std::string& gguP, const std::string& diP);
};