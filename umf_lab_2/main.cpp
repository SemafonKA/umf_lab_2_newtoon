#include <fstream>
#include <iostream>
#include <vector>
#include <array>

#include "ThreeSteppers/Headers/IterSolvers.h"
#include "gaussian_quadrature.h"

// Файл, содержащий в себе пути до файлов, функции f, lambda и gamma
#include "Functions.h"

using namespace std;

#pragma region GLOBAL_OBJECTS
// Глобальная разреженная матрица системы
SparseMatrix global_mat;
// Массив финитных элементов
vector<FiniteElem> finiteElements;
// Массив узлов
vector<Node> nodes;
// Массив сопоставления узлов и первых краевых
vector<S1_node> s1_nodes;
// Массив сопоставления рёбер и вторых краевых
//vector<S23_edge> s2_edges;
//// Массив сопоставления рёбер и третьих краевых
//vector<S23_edge> s3_edges;
#pragma endregion GLOBAL_OBJECTS

void readDataFromFiles() {
   // Считывание данных для структуры узлов nodes
   //auto nodesFile = ifstream(GlobalPaths::nodesPath);
   //if (!nodesFile.is_open())
   //   throw runtime_error("Не удалось открыть файл " + GlobalPaths::nodesPath);
   int size;
   //nodesFile >> size;
   //nodes.resize(size);
   //for (auto& node : nodes) {
   //   nodesFile >> node.x;
   //}
   //nodesFile.close();

   // Считывание данных для структуры прямоугольников finiteElements
   //auto finiteElementsFile = ifstream(GlobalPaths::finiteElementsPath);
   //if (!finiteElementsFile.is_open())
   //   throw runtime_error("Не удалось открыть файл " +
   //      GlobalPaths::finiteElementsPath);
   //finiteElementsFile >> size;
   //finiteElements.resize(size);
   //for (auto& elem : finiteElements) {
   //   finiteElementsFile >> elem.a >> elem.b >> elem.regionNum;
   //}
   //finiteElementsFile.close();

   // Считывание данных для первых краевых условий s1_nodes
   auto s1_nodesFile = ifstream(GlobalPaths::s1_nodesPath);
   if (!s1_nodesFile.is_open())
      throw runtime_error("Не удалось открыть файл " + GlobalPaths::s1_nodesPath);
   s1_nodesFile >> size;
   s1_nodes.resize(size);
   for (auto& s1 : s1_nodes) {
      s1_nodesFile >> s1.node >> s1.funcNum;
   }
   s1_nodesFile.close();

   //// Считывание данных для вторых краевых условий s2_edges
   //auto s2_edgesFile = ifstream(GlobalPaths::s2_edgesPath);
   //if (!s2_edgesFile.is_open())
   //   throw runtime_error("Не удалось открыть файл " + GlobalPaths::s2_edgesPath);
   //s2_edgesFile >> size;
   //s2_edges.resize(size);
   //for (auto& s2 : s2_edges) {
   //   s2_edgesFile >> s2.node1 >> s2.node2 >> s2.funcNum;
   //}
   //s2_edgesFile.close();

   //// Считывание данных для третьих краевых условий s3_edges
   //auto s3_edgesFile = ifstream(GlobalPaths::s3_edgesPath);
   //if (!s3_edgesFile.is_open())
   //   throw runtime_error("Не удалось открыть файл " + GlobalPaths::s3_edgesPath);
   //s3_edgesFile >> size;
   //s3_edges.resize(size);
   //for (auto& s3 : s3_edges) {
   //   s3_edgesFile >> s3.node1 >> s3.node2 >> s3.funcNum;
   //}
   //s3_edgesFile.close();
}

void generatePortrait() {
   global_mat.di.resize(nodes.size());
   global_mat.ig.resize(nodes.size() + 1);

   for (auto& rect : finiteElements) {
      const int elems[2] = { rect.a, rect.b };
      for (int i = 0; i < 2; i++) {
         for (int k = 0; k < i; k++) {
            // Если элемент в верхнем прямоугольнике, то скипаем
            if (elems[k] > elems[i]) continue;

            bool isExist = false;
            // Пробегаем по всей строке для проверки, существует ли такой элемент
            for (auto it = global_mat.ig[elems[i]];
               it < global_mat.ig[elems[i] + 1ll]; it++) {
               if (global_mat.jg[it] == elems[k]) {
                  isExist = true;
                  break;
               }
            }
            if (!isExist) {
               // Ищем, куда вставить элемент портрета
               auto it = global_mat.ig[elems[i]];
               while (it < global_mat.ig[elems[i] + 1ll] &&
                  global_mat.jg[it] < elems[k])
                  it++;

               // Для вставки нужно взять итератор массива от начала, так что...
               global_mat.jg.insert(global_mat.jg.begin() + it, elems[k]);

               // Добавляем всем элементам ig с позиции elems[i]+1 один элемент
               for (auto j = elems[i] + 1; j < global_mat.ig.size(); j++)
                  global_mat.ig[j]++;
            }
         }
      }
   }
   global_mat.ggl.resize(global_mat.jg.size());
   global_mat.ggu.resize(global_mat.jg.size());
}

Matrix getLocalG(const FiniteElem& finiteElem, const std::vector<double>& u, double t) {
   Matrix g = {};

   double x_s = nodes[finiteElem.a].x;
   double hx = abs(nodes[finiteElem.b].x - nodes[finiteElem.a].x);

   double multiplier = lambda_value(finiteElem.regionNum, u[finiteElem.a], nodes[finiteElem.a]);
   multiplier += lambda_value(finiteElem.regionNum, u[finiteElem.b], nodes[finiteElem.b]);
   multiplier = multiplier / 2 / hx;

   constexpr Matrix constMat = {
      1, -1,
      -1, 1
   };
   for (int i = 0; i < finiteSize; i++) {
      for (int j = 0; j < finiteSize; j++) {
         g[i][j] = multiplier * constMat[i][j];
      }
   }

   return g;
}

Matrix getLocalAdditionalMat(const FiniteElem& finiteElem, const std::vector<double>& u, double t) {
   Matrix add = {};
   Matrix dg_1 = {}, dg_2 = {};

   double x_s = nodes[finiteElem.a].x;
   double x_s_1 = nodes[finiteElem.b].x;
   double hx = abs(x_s_1 - x_s);
   double hx_1 = 1.0 / 2.0 / hx;
   double q1 = u[finiteElem.a];
   double q2 = u[finiteElem.b];
   constexpr Matrix constMat = {
       1, -1,
      -1,  1
   };

   // Считаем первую матрицу производных
   double multiplier = lambda_dif_value(finiteElem.regionNum, q1, x_s) / 2 / hx;

   for (int i = 0; i < finiteSize; i++) {
      for (int j = 0; j < finiteSize; j++) {
         dg_1[i][j] = multiplier * constMat[i][j];
      }
   }

   // Считаем вторую матрицу производных
   multiplier = lambda_dif_value(finiteElem.regionNum, q2, x_s_1) / 2 / hx;

   for (int i = 0; i < finiteSize; i++) {
      for (int j = 0; j < finiteSize; j++) {
         dg_2[i][j] = multiplier * constMat[i][j];
      }
   }

   // Считаем добавочную матрицу
   for (int i = 0; i < finiteSize; i++) {
      add[i][0] = dg_1[i][0] * q1 + dg_1[i][1] * q2;
      add[i][1] = dg_2[i][0] * q1 + dg_2[i][1] * q2;
   }

   return add;
}

Matrix getLocalM(const FiniteElem& finiteElem, bool getWithoutSigma = false) {
   Matrix m = {};
   std::function<double(int, double)> maybeSigma;
   if (getWithoutSigma == true) {
      maybeSigma = [](int reg, double x) {
         return 1.0;
      };
   } else {
      maybeSigma = [](int reg, double x) {
         return sigma_value(reg, x);
      };
   }

   double x_s = nodes[finiteElem.a].x;
   double hx = abs(nodes[finiteElem.b].x - nodes[finiteElem.a].x);
   double constCoef = (maybeSigma(finiteElem.regionNum, nodes[finiteElem.a].x) + maybeSigma(finiteElem.regionNum, nodes[finiteElem.b].x)) / 2.0;
   constCoef = constCoef * hx / 6.0;

   static const Matrix constMat = {
      2, 1,
      1, 2
   };

   for (int i = 0; i < finiteSize; i++) {
      for (int j = 0; j < finiteSize; j++) {
         m[i][j] = constCoef * constMat[i][j];
      }
   }

   return m;
}

std::vector<double> getLocalB(const FiniteElem& finiteElem, double t) {
   std::vector<double> b(finiteSize);

   auto M = getLocalM(finiteElem, true);
   for (int i = 0; i < finiteSize; i++) {
      b[i] = 0;
      b[i] += M[i][0] * f_value(finiteElem.regionNum, nodes[finiteElem.a], t);
      b[i] += M[i][1] * f_value(finiteElem.regionNum, nodes[finiteElem.b], t);
   }

   return b;
}

void addLocalMatrixToGlobal(const FiniteElem& finiteElem, SparseMatrix& globalMat, const Matrix& localMat) {
   const int elems[finiteSize] = { finiteElem.a, finiteElem.b };
   for (int i = 0; i < finiteSize; i++) {
      // добавляем все внедиагональные элементы на строке elems[i]
      for (int k = 0; k < i; k++) {
         // Если элемент в верхнем прямоугольнике, то скипаем
         if (elems[k] > elems[i]) {
            continue;
         }

         auto id = globalMat.ig[elems[i]];
         for (id; id < globalMat.ig[elems[i] + 1ll] && globalMat.jg[id] != elems[k]; id++);

         globalMat.ggl[id] += localMat[i][k];
         globalMat.ggu[id] += localMat[k][i];
      }
      // добавляем диагональные элементы
      globalMat.di[elems[i]] += localMat[i][i];
   }
}

void addLocalbToGlobal(const FiniteElem& finiteElem, std::vector<double>& globalVec, const std::vector<double>& localVec) {
   const int elems[finiteSize] = { finiteElem.a, finiteElem.b };
   for (int i = 0; i < finiteSize; i++) {
      globalVec[elems[i]] += localVec[i];
   }
}

//void include_s3() {
//   for (const auto& edge : s3_edges) {
//      double M[2][2] = {};
//      double b[2] = {};
//      double beta = (s3_beta_value(edge.funcNum, nodes[edge.node1]) +
//         s3_beta_value(edge.funcNum, nodes[edge.node2])) /
//         2;
//      double u[2] = { s3_u_value(edge.funcNum, nodes[edge.node1]),
//                     s3_u_value(edge.funcNum, nodes[edge.node2]) };
//      double rp = nodes[edge.node1].r;
//      double hr = nodes[edge.node2].r - nodes[edge.node1].r;
//      double phi_s = nodes[edge.node1].phi;
//      double h_phi = nodes[edge.node2].phi - nodes[edge.node1].phi;
//      // Если краевое задано вдоль оси r
//      if (hr > 1e-7) {
//         M[0][0] = beta * ((hr * hr) / 12 + (rp * hr) / 3);
//         M[0][1] = beta * ((hr * rp) / 6 + (hr * hr) / 12);
//         M[1][0] = M[0][1];
//         M[1][1] = beta * (rp * hr / 3 + hr * hr / 4);
//         b[0] = beta * (u[0] * (hr * rp / 3 + hr * hr / 12) +
//            u[1] * (hr * rp / 6 + hr * hr / 12));
//         b[1] = beta * (u[0] * (hr * rp / 6 + hr * hr / 12) +
//            u[1] * (hr * rp / 3 + hr * hr / 4));
//      }
//      // Если краевое задано вдоль оси phi
//      else {
//         M[0][0] = beta * rp * h_phi / 3;
//         M[0][1] = beta * rp * h_phi / 6;
//         M[1][0] = M[0][1];
//         M[1][1] = beta * rp * h_phi / 3;
//         b[0] = beta * rp * (u[0] * h_phi / 3 + u[1] * h_phi / 6);
//         b[1] = beta * rp * (u[0] * h_phi / 6 + u[1] * h_phi / 3);
//      }
//
//      // добавляем полученный результат в глобальную матрицу
//      int elems[2] = { edge.node1, edge.node2 };
//      for (int i = 0; i < 2; i++) {
//         // добавляем все внедиагональные элементы на строке elems[i]
//         for (int k = 0; k < i; k++) {
//            auto id = global_mat.ig[elems[i]];
//            for (id; id < global_mat.ig[elems[i] + 1ll] &&
//               global_mat.jg[id] != elems[k];
//               id++)
//               ;
//            global_mat.ggl[id] += M[i][k];
//            global_mat.ggu[id] += M[i][k];
//         }
//         // добавляем диагональные элементы и вектор b
//         global_mat.di[elems[i]] += M[i][i];
//         global_b[elems[i]] += b[i];
//      }
//   }
//}
//
//void include_s2() {
//   for (const auto& edge : s2_edges) {
//      double b[2] = {};
//      double theta[2] = { s2_theta_value(edge.funcNum, nodes[edge.node1]),
//                         s2_theta_value(edge.funcNum, nodes[edge.node2]) };
//      double rp = nodes[edge.node1].r;
//      double hr = nodes[edge.node2].r - nodes[edge.node1].r;
//      double phi_s = nodes[edge.node1].phi;
//      double h_phi = nodes[edge.node2].phi - nodes[edge.node1].phi;
//      // Если краевое задано вдоль оси r
//      if (hr > 1e-7) {
//         b[0] = (theta[0] * (hr * rp / 3 + hr * hr / 12) +
//            theta[1] * (hr * rp / 6 + hr * hr / 12));
//         b[1] = (theta[0] * (hr * rp / 6 + hr * hr / 12) +
//            theta[1] * (hr * rp / 3 + hr * hr / 4));
//      }
//      // Если краевое задано вдоль оси phi
//      else {
//         b[0] = rp * (theta[0] * h_phi / 3 + theta[1] * h_phi / 6);
//         b[1] = rp * (theta[0] * h_phi / 6 + theta[1] * h_phi / 3);
//      }
//
//      // добавляем полученный результат в глобальную матрицу
//      int elems[2] = { edge.node1, edge.node2 };
//      for (int i = 0; i < 2; i++) {
//         // добавляем вектор b
//         global_b[elems[i]] += b[i];
//      }
//   }
//}

void include_s1(double t, SparseMatrix& global_mat, vector<double>& global_b) {
   for (const auto& node : s1_nodes) {
      double u = s1_u_value(node.funcNum, nodes[node.node], t);

      // ставим на диагональ значение 1
      global_mat.di[node.node] = 1;
      // ставим в соответствующую ячейку вектора b значение u
      global_b[node.node] = u;
      // зануляем строку в нижнем треугольнике
      for (auto j = global_mat.ig[node.node]; j < global_mat.ig[node.node + 1ll];
         j++) {
         global_mat.ggl[j] = 0;
      }
      // зануляем строку в верхнем треугольнике
      for (int i = node.node + 1; i < global_mat.Size(); i++) {
         for (auto j = global_mat.ig[i]; j < global_mat.ig[i + 1ll]; j++) {
            if (global_mat.jg[j] == node.node) {
               global_mat.ggu[j] = 0;
               break;
            }
         }
      }
   }
}

void getGlobals(
   const SparseMatrix& global_M,
   SparseMatrix& global_mat,
   std::vector<double>& global_d,
   const std::array<std::vector<double>, 2>& slices,
   double t,
   double dt
) {
   std::vector<double> global_b(global_mat.Size());

   // Считаем вектор правой части
   for (auto& elem : finiteElements) {
      auto local_b = getLocalB(elem, t);
      addLocalbToGlobal(elem, global_b, local_b);
   }

   global_d = global_b + (1 / dt) * global_M * slices[0];

   // Считаем матрицу жёсткости с учётом u с начального приближения
   auto global_G = SparseMatrix::copyShape(global_mat);
   for (auto& elem : finiteElements) {
      auto local_g = getLocalG(elem, slices[1], t);
      addLocalMatrixToGlobal(elem, global_G, local_g);
   }

   // Формируем полную глобальную матрицу
   global_mat = global_G + (1 / dt) * global_M;

   // Накладываем краевые
   include_s1(t, global_mat, global_d);
}

void getGlobalsNewtoon(
   const SparseMatrix& global_M,
   SparseMatrix& global_mat,
   std::vector<double>& global_d,
   std::array<std::vector<double>, 2> slices,
   double t,
   double dt,
   double nu
) {
   std::vector<double> global_b(global_mat.Size());

   // Считаем вектор правой части
   for (auto& elem : finiteElements) {
      auto local_b = getLocalB(elem, t);
      addLocalbToGlobal(elem, global_b, local_b);
   }

   global_d = global_b + (1 / dt) * global_M * slices[0];

   // Считаем матрицу жёсткости с учётом u с начального приближения
   auto global_G = SparseMatrix::copyShape(global_mat);
   for (auto& elem : finiteElements) {
      auto local_g = getLocalG(elem, slices[1], t);
      addLocalMatrixToGlobal(elem, global_G, local_g);
   }

   // Формируем полную глобальную матрицу
   global_mat = global_G + (1 / dt) * global_M;

   // Добавочные матрица и вектор Ньютона
   auto global_additional_mat = SparseMatrix::copyShape(global_mat);
   auto global_additional_vec = std::vector<double>(global_mat.Size());
   // Добавляем в глобальную матрицу часть с производными матрицы G
   for (auto& el : finiteElements) {
      auto local_additional_mat = getLocalAdditionalMat(el, slices[1], t);
      addLocalMatrixToGlobal(el, global_additional_mat, local_additional_mat);
   }
   global_additional_vec = global_additional_mat * slices[1];

   global_mat = global_mat + nu * global_additional_mat;
   global_d = global_d + nu * global_additional_vec;

   include_s1(t, global_mat, global_d);
}


int main() {
   setlocale(LC_ALL, "ru-RU");

   // Задаём параметры решателя нелинейного процесса
   constexpr double eps = 1e-13;
   constexpr double maxIter = 500;
   double relaxCoef = 1.0;
   constexpr bool newtoonMode = true;
   if (newtoonMode) {
      cout << "********************\n";
      cout << "** РЕЖИМ НЬЮТОНА  **\n";
      cout << "********************\n\n";
   }

   // Генерируем сетку по пространству
   const double g_begin = 0.0;
   const double g_end = 2.0;
   const double g_step = 1.0 / 4.0;
   nodes.resize(std::floor((g_end - g_begin) / g_step) + 1);
   for (size_t i = 0; i < nodes.size(); i++) {
      nodes[i].x = g_begin + i * g_step;
   }
   finiteElements.resize(nodes.size() - 1);
   for (size_t i = 0; i < nodes.size() - 1; i++) {
      finiteElements[i].a = i;
      finiteElements[i].b = i + 1;
   }
   s1_nodes.resize(2);
   s1_nodes[0].node = 0;
   s1_nodes[1].node = nodes.size() - 1;

   cout << "СЕТКА\n";
   cout << "Начало: " << g_begin << " , конец: " << g_end << ", шаг: " << g_step << endl;
   cout << "Число узлов: " << nodes.size() << endl;
   cout << "Координаты узлов: [" << endl;
   for (auto el : nodes) {
      cout << format("{:20.15e}\n", el.x);
   }
   cout << "]\n\n";

   // Генерируем сетку по времени
   const double t_begin = 0.0;
   const double t_end = 4.0;
   const double t_step = 1.0 / 1.0;
   std::vector<double> t(std::floor((t_end - t_begin) / t_step) + 1);
   for (size_t i = 0; i < t.size(); i++) {
      t[i] = t_begin + i * t_step;
   }

   generatePortrait();

   auto global_M = SparseMatrix::copyShape(global_mat);

   // Считаем глобальную матрицу массы (не будет меняться)
   for (auto& elem : finiteElements) {
      auto local_m = getLocalM(elem, false);
      addLocalMatrixToGlobal(elem, global_M, local_m);
   }

   // Создаём 2 слоя, 1 слой считаем через начальное условие
   std::array<std::vector<double>, 2> slices = {};
   slices[0].resize(global_mat.Size());
   slices[1].resize(global_mat.Size());
   for (size_t i = 0; i < slices[0].size(); i++) {
      slices[1][i] = slices[0][i] = s1_u_value(0, nodes[i], t_begin);
      //slices[1][i] = 0.0;
   }

   // Инициализируем решатель
   IterSolvers::minEps = 1e-15;
   IterSolvers::maxIter = 5000;
   IterSolvers::MSG_Assimetric::Init_LuPrecond(global_mat.Size(), global_mat);
   std::vector<double> errors = { 0.0 };           // Массив погрешностей
   std::vector<double> middle_errors = { 0.0 };    // Массив погрешностей по середине конечных элементов
   std::vector<size_t> number_of_iter = { 0 };     // Массив числа итераций на каждом слое

   for (size_t i = 1; i < t.size(); i++) {
      cout << "Текущее время: " << t[i] << endl;
      double dt = t[i] - t[i - 1];

      std::vector<double> global_d(global_mat.Size());

      // Цикл простой итерации
      double nev = INFINITY;
      size_t j;

      if (!newtoonMode) { // (без Ньютона)
         // Получаем глобальные элементы
         getGlobals(global_M, global_mat, global_d, slices, t[i], dt);
         for (j = 0; j < maxIter && nev > eps; j++) {
            // Решаем СЛАУ
            double ee = 0.0;
            auto newSlice = slices[1]; // Решение q* - временное, копируем его

            IterSolvers::MSG_Assimetric::LuPrecond(global_mat, global_d, newSlice, ee);
            slices[1] = relaxCoef * newSlice + (1.0 - relaxCoef) * (slices[1]);

            // Считаем глобальные матрицы
            getGlobals(global_M, global_mat, global_d, slices, t[i], dt);

            // Находим норму погрешности
            nev = norm(global_mat * slices[1] - global_d) / norm(global_d);
         }
      } else { // С Ньютоном
         for (j = 0; j < maxIter && nev > eps; j++) {
            double nu = 2.0;  // коэффициент дэмпинга
            double newNev = nev;
            auto oldSlice = slices[1];

            // Дэмпинг матрицы на случай, если она перестаёт быть положительно определённой
            while (newNev >= nev && nu > (1.0 / 256.0)) {
               static const vector<double> relaxCoeffs = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };

               nu /= 2.0;
               slices[1] = oldSlice;
               // Получаем глобальные элементы
               getGlobalsNewtoon(global_M, global_mat, global_d, slices, t[i], dt, nu);

               // Решаем СЛАУ
               double ee = 0.0;
               auto newSlice = slices[1]; // Решение q* - временное, копируем его

               IterSolvers::MSG_Assimetric::LuPrecond(global_mat, global_d, newSlice, ee);
               
               // Находим новое приближение функции
               relaxCoef = 1.0;
               slices[1] = relaxCoef * newSlice + (1.0 - relaxCoef) * oldSlice;

               // Считаем глобальные матрицы без добавочных матриц Ньютона
               getGlobals(global_M, global_mat, global_d, slices, t[i], dt);

               // Находим невязку
               newNev = norm(global_mat * slices[1] - global_d) / norm(global_d);

               // Если новая невязка вышла больше, чем была на прошлой итерации, то
               if (newNev > nev) {
                  // Уменьшаем значение коэф. релаксации с шагом 0.1
                  for (int k = relaxCoeffs.size() - 2; k >= 0 && newNev > nev; k--) {
                     slices[1] = relaxCoeffs[k] * newSlice + (1.0 - relaxCoeffs[k]) * oldSlice;
                     // Считаем глобальные матрицы без добавочных матриц Ньютона
                     getGlobals(global_M, global_mat, global_d, slices, t[i], dt);
                     newNev = norm(global_mat * slices[1] - global_d) / norm(global_d);
                  }
               }
            }

            nev = newNev;
         }
      }


      // Выводим чё мы нарешали
      //cout << "Выход из цикла за " << j << " итераций\n";
      number_of_iter.push_back(j);
      cout << "Полученное решение: \n";
      for (auto el : slices[1]) {
         cout << format("{:.15e}\n", el);
      }
      //cout << "Погрешность решения: \n";
      double err = 0.0;
      for (int j = 0; j < slices[1].size(); j++) {
         double tmperr = slices[1][j] - s1_u_value(0, nodes[j], t[i]);
         err += tmperr * tmperr;
         //cout << format("{:5} {:.15e}\n", nodes[j].x, tmperr);
      }
      errors.push_back(std::sqrt(err));

      //cout << "Погрешность решения между узлами: \n";
      err = 0.0;
      for (int j = 0; j < slices[1].size() - 1; j++) {
         double x1 = nodes[j].x;
         double x2 = nodes[j + 1].x;
         double dx = x2 - x1;
         double midx = (x1 + x2) / 2;
         double u = slices[1][j] * (x2 - midx) / dx + slices[1][j + 1] * (midx - x1) / dx;
         double tmperr = u - s1_u_value(0, midx, t[i]);
         err += tmperr * tmperr;
         //cout << format("{:5} {:.15e}\n", midx, tmperr);
      }
      middle_errors.push_back(std::sqrt(err));

      // Сохраняем текущий слой на предыдущий
      slices[0] = slices[1];
      // Если не нужно использовать значение с прошлого слоя как начальное для нового, то
      // раскомментить строчку ниже v
      //slices[1] = std::vector<double>(slices[1].size());
   }

   // Уничтожаем всё, что мы выделили в решателе
   IterSolvers::Destruct();

   cout << "Сводная таблица (по узлам): " << endl;
   cout << "| Время | норма погрешности (узлы) | норма погрешности (середина) | число итераций | \n";
   cout << "| :---: |           :---:          |            :---:             |      :---:     |\n";
   for (int i = 0; i < t.size(); i++) {
      cout << format("| {:^5} | {:^24.15e} | {:^28.15e} | {:^14} |\n", t[i], errors[i], middle_errors[i], number_of_iter[i]);
   }

   return 0;
}
