#include "../Headers/SparseMatrix.h"

using namespace std;

vector<double> ReadVecFromFile(size_t size, const string& path) {
   vector<double> vec(size);
   auto file = ifstream(path);
   if (!file.is_open()) {
      throw runtime_error("Файл " + path + " отсутствует в директории");
   }
   for (size_t i = 0; i < size; i++) {
      file >> vec[i];
   }
   file.close();

   return vec;
}

// Методы матрицы

SparseMatrix SparseMatrix::operator+(const SparseMatrix& other) const {
   SparseMatrix result(*this);
   for (size_t i = 0; i < result.di.size(); i++) {
      result.di.at(i) += other.di.at(i);
   }
   for (size_t i = 0; i < result.ggl.size(); i++) {
      result.ggl.at(i) += other.ggl.at(i);
      result.ggu.at(i) += other.ggu.at(i);
   }
   return result;
}

size_t SparseMatrix::Size() const { return di.size(); }

/// <summary>
/// Умножение матрицы на вектор
/// </summary>
vector<double> SparseMatrix::MultToVec(const vector<double>& right) const {
   vector<double> result(right.size());

   return MultToVec(right, result);
}

/// <summary>
/// Умножение матрицы на вектор
/// </summary>
vector<double>& SparseMatrix::MultToVec(const vector<double>& right, vector<double>& result) const {
   if (right.size() != di.size()) throw runtime_error("Размеры матрицы и вектора не совпадают.");
   if (right.size() != result.size()) throw runtime_error("Размеры матрицы и результирующего вектора не совпадают.");

   for (uint16_t i = 0; i < result.size(); i++) {
      // Умножаем диагональ
      result[i] = di[i] * right[i];

      // Умножаем нижний и верхний треугольники
      for (uint32_t j = ig[i]; j < ig[i + 1ll]; j++) {
         result[i] += ggl[j] * right[jg[j]];
         result[jg[j]] += ggu[j] * right[i];
      }
   }

   return result;
}

vector<double> SparseMatrix::operator*(const vector<double>& right) const { return MultToVec(right); }

SparseMatrix SparseMatrix::multToScalar(double scalar) const {
   SparseMatrix result(*this);

   for (auto& el : result.di) {
      el *= scalar;
   }
   for (size_t i = 0; i < result.ggl.size(); i++) {
      result.ggl[i] *= scalar;
      result.ggu[i] *= scalar;
   }
   return result;
}

/// <summary>
/// Умножение транспонированной матрицы на вектор
/// </summary>
vector<double>& SparseMatrix::TranspMultToVec(const vector<double>& right, vector<double>& result) const {
   if (right.size() != di.size()) throw runtime_error("Размеры матрицы и вектора не совпадают.");
   if (right.size() != result.size()) throw runtime_error("Размеры матрицы и результирующего вектора не совпадают.");

   for (uint16_t i = 0; i < result.size(); i++) {
      // Умножаем диагональ
      result[i] = di[i] * right[i];

      // Умножаем нижний и верхний треугольники
      for (uint32_t j = ig[i]; j < ig[i + 1ll]; j++) {
         result[i] += ggu[j] * right[jg[j]];
         result[jg[j]] += ggl[j] * right[i];
      }
   }

   return result;
}

/// <summary>
/// Умножение транспонированной матрицы на вектор
/// </summary>
vector<double> SparseMatrix::TranspMultToVec(const vector<double>& right) const {
   vector<double> result(right.size());

   return TranspMultToVec(right, result);
}

SparseMatrix& SparseMatrix::operator=(SparseMatrix&& right) noexcept {
   ig = std::move(right.ig);
   jg = std::move(right.jg);
   ggl = std::move(right.ggl);
   ggu = std::move(right.ggu);
   di = std::move(right.di);

   return *this;
}

double SparseMatrix::val(uint16_t row, uint16_t column) {
   // if element lay on diagonal
   if (row == column) return di[row];

   auto& v = row > column ? ggl : ggu;
   if (row < column) std::swap(row, column);

   // find element by his pos
   auto i = ig[row];
   while (i < ig[row + 1ll] && jg[i] < column) i++;

   // if element exists, return him
   if (i < ig[row + 1ll] && jg[i] == column) return v[i];

   // else return zero
   return 0.0;
}

// Конструкторы матрицы

SparseMatrix::SparseMatrix() {}

SparseMatrix::SparseMatrix(const SparseMatrix& right)
    : ig{right.ig}, jg{right.jg}, ggl{right.ggl}, ggu{right.ggu}, di{right.di} {}

// Конструктор перемещения (нужен для метода ReadFromFiles)
SparseMatrix::SparseMatrix(SparseMatrix&& right) noexcept {
   ig = std::move(right.ig);
   jg = std::move(right.jg);
   ggl = std::move(right.ggl);
   ggu = std::move(right.ggu);
   di = std::move(right.di);
}

// Статические методы матрицы

SparseMatrix SparseMatrix::ReadFromFiles(uint16_t matrixSize, const string& igP, const string& jgP, const string& gglP,
                                         const string& gguP, const string& diP) {
   SparseMatrix mat;
   bool isStartFromOne = false;
   {
      mat.ig.resize(matrixSize + 1ll);
      auto igS = ifstream(igP);
      if (!igS.is_open()) throw runtime_error("Файл " + igP + " отсутствует в директории.");
      for (uint16_t i = 0; i <= matrixSize; i++) {
         igS >> mat.ig[i];
      }
      // Если массив ig в файле начинался с 1, то меняем его под наши параметры (под 0)
      if (isStartFromOne = mat.ig[0]) {
         for (uint16_t i = 0; i <= matrixSize; i++) {
            mat.ig[i]--;
         }
      }
   }

   {
      auto jgS = ifstream(jgP);
      if (!jgS.is_open()) throw runtime_error("Файл " + jgP + " отсутствует в директории.");
      mat.jg.resize(mat.ig.back());
      for (uint32_t i = 0; i < mat.jg.size(); i++) {
         jgS >> mat.jg[i];
         if (isStartFromOne) {
            mat.jg[i]--;
         }
      }
   }
   try {
      mat.di = ReadVecFromFile(matrixSize, diP);
      mat.ggl = ReadVecFromFile(mat.jg.size(), gglP);
      mat.ggu = ReadVecFromFile(mat.jg.size(), gguP);
   } catch (exception& e) {
      throw e;
   }

   return mat;
}

SparseMatrix SparseMatrix::copyShape(const SparseMatrix& other) {
   SparseMatrix mat;
   mat.ig = other.ig;
   mat.jg = other.jg;
   mat.di.resize(other.di.size());
   mat.ggl.resize(other.ggl.size());
   mat.ggu.resize(other.ggu.size());
   return mat;
}