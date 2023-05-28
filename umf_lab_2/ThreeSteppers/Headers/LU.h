#pragma once
#include <vector>
#include "SparseMatrix.h"

// Неполное разложение LU(sq) матрицы разреженного строчно-столбцового формата SparseMatrix
// Не хранит портрет матрицы, но использует портрет исходной матрицы (а также ссылается на неё)
class LU {
// Блок внутренних переменных разложения LU
public:
   const SparseMatrix* parent = nullptr;

   // Вектор диагональных элементов LU разложения. В данном случае диагонали L и U совпадают
   std::vector<double> di;

   // Вектор элементов нижнего треугольника L
   std::vector<double> ggl;

   // Вектор элементов верхнего треугольника U
   std::vector<double> ggu;

// Блок основных конструкторов класса
public:

   /// <summary>
   /// Конструктор с резервированием памяти под разложение
   /// </summary>
   /// <param name="diSize"> - размер диагонали,</param>
   /// <param name="luSize"> - размер массивов нижнего и верхнего треугольника</param>
   LU(size_t diSize, size_t luSize);

   /// <summary>
   /// Конструктор с построением неполного LU(sq)-разложения по матрице mat
   /// </summary>
   /// <param name="mat"> - матрица, по которой построится LU-разложение, с привязкой этой матрицы к объекту</param>
   LU(const SparseMatrix& mat);

// Блок основных нестатических методов класса
public:

   /// <summary>
   /// Разложить матрицу mat в неполное LU(sq) - разложение
   /// </summary>
   /// <param name="mat"> - матрица, которую требуется разложить. Она же будет использоваться для просмотра портрета матриц</param>
   void MakeLuFor(const SparseMatrix& mat);

   /// <summary>
   /// Метод изменения размера разложения
   /// </summary>
   /// <param name="diSize"> - размер диагонали,</param>
   /// <param name="luSize"> - размер массивов нижнего и верхнего треугольника</param>
   void Resize(size_t diSize, size_t luSize);

// Умножение матриц на вектор

   /// <summary>
   /// Умножение нижней матрицы L на вектор vec. Выделяет память под вектор ответа, не меняет матрицу LU
   /// </summary>
   /// <param name="vec"> - вектор, на который будет происходить умножение матрицы L</param>
   /// <returns>Вектор с результатом перемножения (выделяется в памяти)</returns>
   std::vector<double> LMultToVec(const std::vector<double>& vec) const;

   /// <summary>
   /// Умножение нижней матрицы L на вектор vec. Ответ записывается в вектор ans, не меняет матрицу LU
   /// </summary>
   /// <param name="vec"> - вектор, на который будет происходить умножение матрицы L;</param>
   /// <param name="ans"> - вектор, куда запишется ответ без выделения памяти (должен отличаться от vec!)</param>
   /// <returns>Ссылка на вектор ans</returns>
   std::vector<double>& LMultToVec(const std::vector<double>& vec, std::vector<double>& ans) const;

   /// <summary>
   /// Умножение нижней матрицы L^T на вектор vec. Выделяет память под вектор ответа, не меняет матрицу LU
   /// </summary>
   /// <param name="vec"> - вектор, на который будет происходить умножение матрицы L</param>
   /// <returns>Вектор с результатом перемножения (выделяется в памяти)</returns>
   std::vector<double> LTranspMultToVec(const std::vector<double>& vec) const;

   /// <summary>
   /// Умножение нижней матрицы L^T на вектор vec. Ответ записывается в вектор ans, не меняет матрицу LU
   /// </summary>
   /// <param name="vec"> - вектор, на который будет происходить умножение матрицы L;</param>
   /// <param name="ans"> - вектор, куда запишется ответ без выделения памяти (должен отличаться от vec!)</param>
   /// <returns>Ссылка на вектор ans</returns>
   std::vector<double>& LTranspMultToVec(const std::vector<double>& vec, std::vector<double>& ans) const;

   /// <summary>
   /// Умножение верхней матрицы U на вектор vec. Выделяет память под вектор ответа, не меняет матрицу LU
   /// </summary>
   /// <param name="vec"> - вектор, на который будет происходить умножение матрицы U</param>
   /// <returns>Вектор с результатом перемножения (выделяется в памяти)</returns>
   std::vector<double> UMultToVec(const std::vector<double>& vec) const;

   /// <summary>
   /// Умножение верхней матрицы U на вектор vec. Ответ записывается в вектор ans, не меняет матрицу LU
   /// </summary>
   /// <param name="vec"> - вектор, на который будет происходить умножение матрицы U;</param>
   /// <param name="ans"> - вектор, куда запишется ответ без выделения памяти (должен отличаться от vec!)</param>
   /// <returns>Ссылка на вектор ans</returns>
   std::vector<double>& UMultToVec(const std::vector<double>& vec, std::vector<double>& ans) const;

   /// <summary>
   /// Умножение верхней матрицы U^T на вектор vec. Выделяет память под вектор ответа, не меняет матрицу LU
   /// </summary>
   /// <param name="vec"> - вектор, на который будет происходить умножение матрицы U</param>
   /// <returns>Вектор с результатом перемножения (выделяется в памяти)</returns>
   std::vector<double> UTranspMultToVec(const std::vector<double>& vec) const;

   /// <summary>
   /// Умножение верхней матрицы U^T на вектор vec. Ответ записывается в вектор ans, не меняет матрицу LU
   /// </summary>
   /// <param name="vec"> - вектор, на который будет происходить умножение матрицы U;</param>
   /// <param name="ans"> - вектор, куда запишется ответ без выделения памяти (должен отличаться от vec!)</param>
   /// <returns>Ссылка на вектор ans</returns>
   std::vector<double>& UTranspMultToVec(const std::vector<double>& vec, std::vector<double>& ans) const;

// Решение слау с использованием матриц и вектора правой части

   /// <summary>
   /// Решение слау вида Lx = right. Не выделяет память под вектор x, не меняет матрицы LU
   /// </summary>
   /// <param name="right"> - вектор правой части уравнения;</param>
   /// <param name="x"> - вектор, куда будет записан ответ. Должен быть с уже выделенной памятью. Должен отличаться от right!</param>
   /// <returns>ссылка на вектор x</returns>
   std::vector<double>& LSlauSolve(const std::vector<double>& right, std::vector<double>& x) const;

   /// <summary>
   /// Решение слау вида Lx = right. Выделяет память под вектор x, не меняет матрицы LU
   /// </summary>
   /// <param name="right"> - вектор правой части уравнения;</param>
   /// <returns>Полученный вектор x</returns>
   std::vector<double> LSlauSolve(const std::vector<double>& right) const;

   /// <summary>
   /// Решение слау вида L^T * x = right. Не выделяет память под вектор x, не меняет матрицы LU
   /// </summary>
   /// <param name="right"> - вектор правой части уравнения;</param>
   /// <param name="x"> - вектор, куда будет записан ответ. Должен быть с уже выделенной памятью. Должен отличаться от right!</param>
   /// <returns>ссылка на вектор x</returns>
   std::vector<double>& LTranspSlauSolve(const std::vector<double>& right, std::vector<double>& x) const;

   /// <summary>
   /// Решение слау вида L^T * x = right. Выделяет память под вектор x, не меняет матрицы LU
   /// </summary>
   /// <param name="right"> - вектор правой части уравнения;</param>
   /// <returns>Полученный вектор x</returns>
   std::vector<double> LTranspSlauSolve(const std::vector<double>& right) const;

   /// <summary>
   /// Решение слау вида Ux = right. Не выделяет память под вектор x, не меняет матрицы LU
   /// </summary>
   /// <param name="right"> - вектор правой части уравнения;</param>
   /// <param name="x"> - вектор, куда будет записан ответ. Должен быть с уже выделенной памятью. Должен отличаться от right!</param>
   /// <returns>ссылка на вектор x</returns>
   std::vector<double>& USlauSolve(const std::vector<double>& right, std::vector<double>& x) const;

   /// <summary>
   /// Решение слау вида Ux = right. Выделяет память под вектор x, не меняет матрицы LU
   /// </summary>
   /// <param name="right"> - вектор правой части уравнения;</param>
   /// <returns>Полученный вектор x</returns>
   std::vector<double> USlauSolve(const std::vector<double>& right) const;

   /// <summary>
   /// Решение слау вида U^T * x = right. Не выделяет память под вектор x, не меняет матрицы LU
   /// </summary>
   /// <param name="right"> - вектор правой части уравнения;</param>
   /// <param name="x"> - вектор, куда будет записан ответ. Должен быть с уже выделенной памятью. Должен отличаться от right!</param>
   /// <returns>ссылка на вектор x</returns>
   std::vector<double>& UTranspSlauSolve(const std::vector<double>& right, std::vector<double>& x) const;

   /// <summary>
   /// Решение слау вида U^T * x = right. Выделяет память под вектор x, не меняет матрицы LU
   /// </summary>
   /// <param name="right"> - вектор правой части уравнения;</param>
   /// <returns>Полученный вектор x</returns>
   std::vector<double> UTranspSlauSolve(const std::vector<double>& right) const;
};