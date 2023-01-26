#include <vector>
#include <random>
#include <DistanceMatrixDistorsion.h>

ttk::DistanceMatrixDistorsion::DistanceMatrixDistorsion() {
  this->setDebugMsgPrefix("DistanceMatrixDistorsion");
}


int ttk::DistanceMatrixDistorsion::test(int n, std::vector<double> &res) const
{

  double sim;
  std::vector<std::vector<double>> mat1 =
{
{0, 3, 30, 11, 20, 21, 37, 22, 35, 16, 28, 35, 43, 14, 30, 18, 44, 49, 30, 46, 30, 10, 28, 44, 14, 11, 18, 13, 13, 30, 14, 4, 3, 33, 20, 7, 48, 1, 27, 17, 50, 18, 9, 11, 17, 2, 14, 47, 42, 34},
{3, 0, 49, 11, 19, 16, 38, 2, 11, 48, 20, 31, 28, 3, 5, 35, 2, 13, 3, 47, 24, 22, 5, 35, 33, 2, 20, 23, 7, 11, 50, 37, 50, 17, 4, 22, 37, 32, 35, 12, 34, 0, 4, 15, 17, 16, 41, 40, 31, 13},
{30, 49, 0, 27, 24, 43, 15, 1, 32, 31, 33, 38, 27, 35, 11, 25, 13, 23, 0, 16, 31, 2, 33, 17, 0, 7, 46, 9, 44, 45, 8, 13, 17, 23, 30, 41, 34, 20, 16, 5, 13, 44, 47, 30, 39, 15, 24, 14, 18, 36},
{11, 11, 27, 0, 49, 31, 10, 43, 13, 7, 3, 42, 21, 10, 49, 42, 9, 20, 16, 11, 37, 42, 33, 38, 1, 46, 30, 40, 30, 33, 34, 9, 41, 0, 18, 4, 38, 29, 1, 16, 13, 9, 39, 33, 42, 39, 33, 13, 15, 8},
{20, 19, 24, 49, 0, 30, 7, 12, 45, 1, 21, 26, 23, 16, 50, 19, 33, 24, 49, 30, 0, 24, 15, 25, 39, 43, 50, 44, 16, 7, 28, 22, 3, 14, 38, 14, 13, 6, 39, 46, 44, 29, 30, 30, 28, 5, 35, 25, 45, 19},
{21, 16, 43, 31, 30, 0, 23, 21, 34, 23, 43, 45, 26, 42, 11, 5, 26, 50, 1, 45, 43, 29, 24, 1, 9, 28, 8, 40, 14, 14, 6, 0, 21, 50, 35, 13, 44, 47, 2, 34, 28, 37, 32, 40, 10, 17, 22, 42, 24, 16},
{37, 38, 15, 10, 7, 23, 0, 24, 39, 50, 31, 3, 31, 0, 11, 14, 29, 23, 50, 12, 20, 19, 47, 37, 33, 1, 5, 6, 37, 9, 11, 43, 8, 8, 7, 6, 12, 15, 23, 1, 5, 15, 9, 1, 12, 31, 43, 39, 31, 16},
{22, 2, 1, 43, 12, 21, 24, 0, 15, 45, 9, 39, 9, 23, 0, 3, 39, 13, 40, 32, 45, 13, 11, 4, 35, 16, 43, 39, 5, 33, 35, 12, 41, 18, 12, 44, 4, 45, 13, 43, 5, 20, 16, 14, 17, 9, 30, 28, 19, 8},
{35, 11, 32, 13, 45, 34, 39, 15, 0, 41, 16, 8, 30, 28, 35, 28, 12, 10, 33, 49, 46, 16, 10, 2, 27, 38, 12, 49, 25, 26, 9, 42, 17, 35, 19, 30, 23, 24, 25, 34, 27, 15, 8, 21, 37, 32, 12, 20, 9, 47},
{16, 48, 31, 7, 1, 23, 50, 45, 41, 0, 15, 46, 49, 3, 6, 1, 34, 30, 43, 31, 12, 40, 16, 45, 28, 47, 2, 48, 22, 15, 4, 4, 0, 4, 37, 29, 43, 28, 6, 7, 12, 22, 11, 20, 11, 3, 17, 30, 28, 2},
{28, 20, 33, 3, 21, 43, 31, 9, 16, 15, 0, 46, 7, 1, 17, 32, 49, 40, 44, 17, 18, 36, 9, 36, 32, 32, 46, 14, 35, 44, 7, 16, 22, 4, 30, 50, 29, 49, 14, 49, 35, 18, 28, 9, 16, 40, 41, 32, 1, 46},
{35, 31, 38, 42, 26, 45, 3, 39, 8, 46, 46, 0, 50, 43, 46, 0, 20, 42, 22, 31, 26, 30, 18, 9, 44, 22, 34, 47, 14, 16, 48, 11, 16, 2, 43, 18, 49, 7, 5, 12, 50, 38, 40, 38, 10, 18, 23, 11, 35, 20},
{43, 28, 27, 21, 23, 26, 31, 9, 30, 49, 7, 50, 0, 14, 33, 22, 43, 37, 5, 44, 33, 5, 30, 8, 46, 33, 1, 34, 0, 0, 20, 6, 42, 38, 26, 4, 31, 50, 48, 48, 16, 29, 49, 21, 27, 0, 38, 3, 12, 38},
{14, 3, 35, 10, 16, 42, 0, 23, 28, 3, 1, 43, 14, 0, 29, 11, 22, 38, 29, 25, 37, 30, 47, 12, 16, 49, 2, 16, 38, 29, 31, 19, 47, 34, 10, 3, 16, 21, 44, 32, 12, 10, 17, 26, 16, 17, 20, 38, 28, 26},
{30, 5, 11, 49, 50, 11, 11, 0, 35, 6, 17, 46, 33, 29, 0, 32, 30, 24, 27, 33, 30, 31, 32, 35, 35, 23, 25, 9, 14, 15, 44, 32, 39, 5, 16, 28, 19, 43, 19, 34, 40, 46, 23, 19, 12, 16, 47, 20, 31, 2},
{18, 35, 25, 42, 19, 5, 14, 3, 28, 1, 32, 0, 22, 11, 32, 0, 20, 0, 16, 29, 17, 20, 38, 46, 39, 47, 22, 1, 25, 45, 21, 19, 49, 30, 5, 22, 11, 26, 4, 49, 43, 34, 25, 1, 25, 40, 12, 14, 13, 42},
{44, 2, 13, 9, 33, 26, 29, 39, 12, 34, 49, 20, 43, 22, 30, 20, 0, 30, 22, 46, 25, 36, 49, 34, 5, 30, 40, 28, 8, 6, 3, 37, 40, 45, 16, 41, 21, 18, 22, 34, 5, 7, 4, 21, 40, 47, 7, 28, 37, 4},
{49, 13, 23, 20, 24, 50, 23, 13, 10, 30, 40, 42, 37, 38, 24, 0, 30, 0, 4, 40, 12, 18, 2, 17, 32, 30, 13, 26, 25, 3, 14, 26, 49, 7, 16, 21, 10, 28, 26, 46, 14, 1, 1, 20, 35, 20, 41, 36, 49, 25},
{30, 3, 0, 16, 49, 1, 50, 40, 33, 43, 44, 22, 5, 29, 27, 16, 22, 4, 0, 41, 0, 12, 13, 43, 27, 37, 18, 3, 46, 23, 49, 24, 27, 19, 43, 50, 22, 29, 11, 46, 21, 41, 44, 22, 22, 21, 21, 0, 3, 11},
{46, 47, 16, 11, 30, 45, 12, 32, 49, 31, 17, 31, 44, 25, 33, 29, 46, 40, 41, 0, 15, 7, 29, 9, 9, 43, 10, 20, 29, 13, 44, 47, 24, 43, 5, 31, 27, 41, 40, 23, 13, 9, 23, 41, 13, 38, 41, 18, 21, 36},
{30, 24, 31, 37, 0, 43, 20, 45, 46, 12, 18, 26, 33, 37, 30, 17, 25, 12, 0, 15, 0, 31, 15, 5, 44, 21, 42, 41, 17, 8, 21, 5, 29, 32, 6, 22, 46, 7, 19, 45, 9, 2, 31, 39, 37, 49, 40, 40, 21, 15},
{10, 22, 2, 42, 24, 29, 19, 13, 16, 40, 36, 30, 5, 30, 31, 20, 36, 18, 12, 7, 31, 0, 8, 6, 27, 35, 40, 23, 47, 45, 46, 2, 42, 43, 36, 43, 15, 39, 41, 2, 50, 3, 3, 41, 29, 27, 36, 10, 43, 6},
{28, 5, 33, 33, 15, 24, 47, 11, 10, 16, 9, 18, 30, 47, 32, 38, 49, 2, 13, 29, 15, 8, 0, 18, 27, 2, 25, 19, 38, 8, 25, 46, 42, 10, 1, 39, 1, 44, 35, 50, 38, 50, 1, 39, 33, 43, 28, 45, 22, 23},
{44, 35, 17, 38, 25, 1, 37, 4, 2, 45, 36, 9, 8, 12, 35, 46, 34, 17, 43, 9, 5, 6, 18, 0, 1, 27, 50, 5, 22, 46, 20, 47, 24, 22, 14, 30, 12, 18, 10, 32, 12, 43, 33, 38, 33, 28, 8, 25, 14, 0},
{14, 33, 0, 1, 39, 9, 33, 35, 27, 28, 32, 44, 46, 16, 35, 39, 5, 32, 27, 9, 44, 27, 27, 1, 0, 9, 31, 3, 40, 46, 33, 26, 32, 0, 43, 46, 44, 33, 7, 20, 44, 37, 43, 44, 19, 44, 0, 33, 38, 40},
{11, 2, 7, 46, 43, 28, 1, 16, 38, 47, 32, 22, 33, 49, 23, 47, 30, 30, 37, 43, 21, 35, 2, 27, 9, 0, 30, 7, 21, 28, 49, 3, 18, 50, 16, 2, 30, 23, 27, 4, 3, 8, 21, 12, 23, 49, 28, 23, 40, 25},
{18, 20, 46, 30, 50, 8, 5, 43, 12, 2, 46, 34, 1, 2, 25, 22, 40, 13, 18, 10, 42, 40, 25, 50, 31, 30, 0, 2, 50, 20, 3, 17, 4, 10, 10, 45, 46, 15, 41, 5, 23, 19, 2, 29, 28, 26, 28, 42, 18, 14},
{13, 23, 9, 40, 44, 40, 6, 39, 49, 48, 14, 47, 34, 16, 9, 1, 28, 26, 3, 20, 41, 23, 19, 5, 3, 7, 2, 0, 34, 33, 38, 27, 50, 16, 35, 46, 46, 4, 8, 32, 33, 35, 13, 27, 38, 18, 31, 34, 7, 41},
{13, 7, 44, 30, 16, 14, 37, 5, 25, 22, 35, 14, 0, 38, 14, 25, 8, 25, 46, 29, 17, 47, 38, 22, 40, 21, 50, 34, 0, 31, 45, 17, 49, 16, 27, 18, 35, 42, 46, 41, 2, 18, 3, 39, 26, 8, 23, 35, 10, 24},
{30, 11, 45, 33, 7, 14, 9, 33, 26, 15, 44, 16, 0, 29, 15, 45, 6, 3, 23, 13, 8, 45, 8, 46, 46, 28, 20, 33, 31, 0, 48, 24, 28, 13, 1, 10, 15, 4, 26, 30, 2, 36, 13, 3, 0, 49, 4, 10, 50, 13},
{14, 50, 8, 34, 28, 6, 11, 35, 9, 4, 7, 48, 20, 31, 44, 21, 3, 14, 49, 44, 21, 46, 25, 20, 33, 49, 3, 38, 45, 48, 0, 15, 47, 48, 44, 8, 4, 46, 15, 7, 26, 38, 47, 29, 27, 7, 41, 1, 18, 33},
{4, 37, 13, 9, 22, 0, 43, 12, 42, 4, 16, 11, 6, 19, 32, 19, 37, 26, 24, 47, 5, 2, 46, 47, 26, 3, 17, 27, 17, 24, 15, 0, 27, 37, 29, 30, 22, 11, 3, 22, 49, 19, 32, 50, 9, 37, 21, 5, 7, 29},
{3, 50, 17, 41, 3, 21, 8, 41, 17, 0, 22, 16, 42, 47, 39, 49, 40, 49, 27, 24, 29, 42, 42, 24, 32, 18, 4, 50, 49, 28, 47, 27, 0, 29, 9, 28, 36, 40, 40, 3, 4, 44, 44, 33, 9, 34, 28, 0, 27, 42},
{33, 17, 23, 0, 14, 50, 8, 18, 35, 4, 4, 2, 38, 34, 5, 30, 45, 7, 19, 43, 32, 43, 10, 22, 0, 50, 10, 16, 16, 13, 48, 37, 29, 0, 50, 25, 36, 47, 9, 47, 32, 30, 22, 38, 8, 42, 41, 15, 50, 3},
{20, 4, 30, 18, 38, 35, 7, 12, 19, 37, 30, 43, 26, 10, 16, 5, 16, 16, 43, 5, 6, 36, 1, 14, 43, 16, 10, 35, 27, 1, 44, 29, 9, 50, 0, 43, 17, 5, 43, 18, 19, 36, 34, 50, 20, 5, 0, 17, 1, 3},
{7, 22, 41, 4, 14, 13, 6, 44, 30, 29, 50, 18, 4, 3, 28, 22, 41, 21, 50, 31, 22, 43, 39, 30, 46, 2, 45, 46, 18, 10, 8, 30, 28, 25, 43, 0, 27, 18, 38, 45, 27, 44, 6, 33, 7, 23, 19, 7, 6, 0},
{48, 37, 34, 38, 13, 44, 12, 4, 23, 43, 29, 49, 31, 16, 19, 11, 21, 10, 22, 27, 46, 15, 1, 12, 44, 30, 46, 46, 35, 15, 4, 22, 36, 36, 17, 27, 0, 3, 49, 4, 36, 40, 12, 45, 41, 27, 37, 48, 2, 43},
{1, 32, 20, 29, 6, 47, 15, 45, 24, 28, 49, 7, 50, 21, 43, 26, 18, 28, 29, 41, 7, 39, 44, 18, 33, 23, 15, 4, 42, 4, 46, 11, 40, 47, 5, 18, 3, 0, 17, 6, 36, 6, 4, 15, 15, 3, 9, 16, 14, 12},
{27, 35, 16, 1, 39, 2, 23, 13, 25, 6, 14, 5, 48, 44, 19, 4, 22, 26, 11, 40, 19, 41, 35, 10, 7, 27, 41, 8, 46, 26, 15, 3, 40, 9, 43, 38, 49, 17, 0, 17, 43, 29, 39, 13, 32, 33, 39, 20, 41, 23},
{17, 12, 5, 16, 46, 34, 1, 43, 34, 7, 49, 12, 48, 32, 34, 49, 34, 46, 46, 23, 45, 2, 50, 32, 20, 4, 5, 32, 41, 30, 7, 22, 3, 47, 18, 45, 4, 6, 17, 0, 24, 40, 10, 41, 5, 1, 20, 1, 38, 12},
{50, 34, 13, 13, 44, 28, 5, 5, 27, 12, 35, 50, 16, 12, 40, 43, 5, 14, 21, 13, 9, 50, 38, 12, 44, 3, 23, 33, 2, 2, 26, 49, 4, 32, 19, 27, 36, 36, 43, 24, 0, 8, 48, 46, 30, 34, 24, 16, 7, 22},
{18, 0, 44, 9, 29, 37, 15, 20, 15, 22, 18, 38, 29, 10, 46, 34, 7, 1, 41, 9, 2, 3, 50, 43, 37, 8, 19, 35, 18, 36, 38, 19, 44, 30, 36, 44, 40, 6, 29, 40, 8, 0, 38, 23, 47, 33, 14, 30, 11, 6},
{9, 4, 47, 39, 30, 32, 9, 16, 8, 11, 28, 40, 49, 17, 23, 25, 4, 1, 44, 23, 31, 3, 1, 33, 43, 21, 2, 13, 3, 13, 47, 32, 44, 22, 34, 6, 12, 4, 39, 10, 48, 38, 0, 8, 32, 6, 22, 50, 22, 3},
{11, 15, 30, 33, 30, 40, 1, 14, 21, 20, 9, 38, 21, 26, 19, 1, 21, 20, 22, 41, 39, 41, 39, 38, 44, 12, 29, 27, 39, 3, 29, 50, 33, 38, 50, 33, 45, 15, 13, 41, 46, 23, 8, 0, 44, 14, 33, 6, 15, 17},
{17, 17, 39, 42, 28, 10, 12, 17, 37, 11, 16, 10, 27, 16, 12, 25, 40, 35, 22, 13, 37, 29, 33, 33, 19, 23, 28, 38, 26, 0, 27, 9, 9, 8, 20, 7, 41, 15, 32, 5, 30, 47, 32, 44, 0, 21, 1, 47, 5, 42},
{2, 16, 15, 39, 5, 17, 31, 9, 32, 3, 40, 18, 0, 17, 16, 40, 47, 20, 21, 38, 49, 27, 43, 28, 44, 49, 26, 18, 8, 49, 7, 37, 34, 42, 5, 23, 27, 3, 33, 1, 34, 33, 6, 14, 21, 0, 16, 10, 49, 5},
{14, 41, 24, 33, 35, 22, 43, 30, 12, 17, 41, 23, 38, 20, 47, 12, 7, 41, 21, 41, 40, 36, 28, 8, 0, 28, 28, 31, 23, 4, 41, 21, 28, 41, 0, 19, 37, 9, 39, 20, 24, 14, 22, 33, 1, 16, 0, 36, 48, 16},
{47, 40, 14, 13, 25, 42, 39, 28, 20, 30, 32, 11, 3, 38, 20, 14, 28, 36, 0, 18, 40, 10, 45, 25, 33, 23, 42, 34, 35, 10, 1, 5, 0, 15, 17, 7, 48, 16, 20, 1, 16, 30, 50, 6, 47, 10, 36, 0, 14, 20},
{42, 31, 18, 15, 45, 24, 31, 19, 9, 28, 1, 35, 12, 28, 31, 13, 37, 49, 3, 21, 21, 43, 22, 14, 38, 40, 18, 7, 10, 50, 18, 7, 27, 50, 1, 6, 2, 14, 41, 38, 7, 11, 22, 15, 5, 49, 48, 14, 0, 47},
{34, 13, 36, 8, 19, 16, 16, 8, 47, 2, 46, 20, 38, 26, 2, 42, 4, 25, 11, 36, 15, 6, 23, 0, 40, 25, 14, 41, 24, 13, 33, 29, 42, 3, 3, 0, 43, 12, 23, 12, 22, 6, 3, 17, 42, 5, 16, 20, 47, 0}};
 std::vector<std::vector<double>> mat2 =
{
{0, 2, 12, 50, 13, 4, 44, 36, 36, 25, 34, 40, 22, 33, 21, 29, 16, 40, 39, 33, 35, 9, 0, 26, 3, 18, 17, 41, 4, 19, 49, 35, 19, 24, 24, 40, 25, 16, 19, 50, 5, 5, 5, 49, 7, 26, 43, 21, 9, 43},
{2, 0, 43, 18, 15, 35, 19, 29, 44, 4, 49, 39, 7, 43, 10, 4, 41, 23, 28, 39, 35, 35, 47, 28, 0, 35, 17, 49, 19, 34, 11, 15, 12, 27, 6, 23, 13, 11, 46, 6, 22, 13, 28, 31, 2, 41, 41, 20, 24, 16},
{12, 43, 0, 39, 15, 20, 34, 47, 10, 25, 43, 50, 41, 40, 16, 50, 6, 20, 12, 17, 21, 33, 45, 23, 7, 3, 43, 26, 36, 27, 0, 25, 4, 46, 10, 28, 15, 11, 11, 3, 16, 3, 36, 7, 24, 13, 8, 35, 47, 41},
{50, 18, 39, 0, 50, 30, 1, 29, 24, 50, 16, 12, 11, 29, 12, 17, 14, 35, 48, 47, 41, 11, 8, 14, 47, 17, 15, 2, 26, 27, 42, 9, 26, 46, 41, 46, 9, 10, 36, 1, 17, 10, 49, 46, 2, 7, 7, 26, 5, 31},
{13, 15, 15, 50, 0, 8, 38, 29, 28, 34, 7, 8, 3, 25, 1, 28, 3, 18, 49, 46, 42, 26, 36, 22, 46, 32, 20, 35, 18, 29, 8, 41, 17, 7, 1, 37, 13, 18, 26, 0, 39, 29, 15, 14, 17, 37, 27, 47, 32, 30},
{4, 35, 20, 30, 8, 0, 23, 5, 46, 7, 21, 18, 3, 27, 30, 15, 13, 20, 16, 5, 22, 14, 45, 29, 16, 28, 22, 30, 29, 15, 37, 42, 30, 29, 1, 33, 43, 3, 20, 2, 49, 39, 50, 47, 10, 45, 15, 25, 45, 8},
{44, 19, 34, 1, 38, 23, 0, 30, 26, 48, 2, 5, 26, 2, 33, 11, 21, 26, 47, 15, 46, 21, 37, 15, 0, 6, 5, 5, 23, 2, 49, 1, 45, 43, 18, 44, 24, 29, 40, 24, 34, 17, 12, 37, 38, 22, 7, 17, 49, 34},
{36, 29, 47, 29, 29, 5, 30, 0, 35, 26, 8, 18, 42, 47, 34, 38, 49, 40, 32, 35, 40, 38, 40, 47, 26, 24, 12, 39, 30, 40, 13, 7, 6, 16, 27, 13, 1, 50, 6, 7, 27, 22, 31, 50, 44, 32, 44, 12, 25, 23},
{36, 44, 10, 24, 28, 46, 26, 35, 0, 12, 0, 41, 36, 15, 44, 33, 34, 44, 19, 18, 18, 26, 20, 15, 16, 2, 6, 48, 7, 4, 20, 15, 16, 3, 46, 9, 22, 42, 49, 14, 46, 15, 36, 30, 44, 30, 23, 40, 7, 13},
{25, 4, 25, 50, 34, 7, 48, 26, 12, 0, 0, 25, 23, 22, 9, 30, 10, 32, 36, 23, 10, 37, 27, 47, 47, 17, 38, 14, 38, 14, 39, 4, 12, 17, 43, 18, 22, 3, 50, 18, 38, 8, 37, 1, 33, 22, 30, 9, 14, 32},
{34, 49, 43, 16, 7, 21, 2, 8, 0, 0, 0, 20, 49, 44, 13, 21, 30, 48, 2, 37, 2, 21, 24, 16, 1, 26, 47, 47, 5, 1, 0, 39, 21, 16, 13, 8, 22, 24, 50, 2, 35, 5, 21, 9, 19, 27, 1, 1, 45, 3},
{40, 39, 50, 12, 8, 18, 5, 18, 41, 25, 20, 0, 40, 42, 31, 26, 39, 44, 28, 22, 45, 8, 42, 2, 33, 5, 27, 11, 29, 15, 3, 5, 5, 49, 44, 47, 19, 36, 1, 16, 26, 9, 21, 44, 13, 29, 2, 41, 44, 33},
{22, 7, 41, 11, 3, 3, 26, 42, 36, 23, 49, 40, 0, 18, 26, 0, 47, 9, 5, 50, 23, 25, 46, 6, 33, 42, 7, 7, 48, 23, 38, 46, 3, 41, 15, 42, 28, 15, 18, 40, 45, 3, 1, 43, 36, 24, 17, 8, 21, 30},
{33, 43, 40, 29, 25, 27, 2, 47, 15, 22, 44, 42, 18, 0, 27, 31, 17, 38, 42, 35, 17, 5, 34, 43, 32, 46, 5, 37, 47, 39, 40, 10, 46, 46, 1, 45, 41, 34, 23, 50, 27, 27, 19, 46, 1, 30, 46, 14, 30, 14},
{21, 10, 16, 12, 1, 30, 33, 34, 44, 9, 13, 31, 26, 27, 0, 19, 38, 14, 19, 37, 13, 37, 5, 30, 21, 2, 32, 49, 43, 38, 23, 35, 8, 3, 30, 40, 2, 4, 14, 0, 1, 40, 31, 41, 22, 0, 9, 3, 0, 37},
{29, 4, 50, 17, 28, 15, 11, 38, 33, 30, 21, 26, 0, 31, 19, 0, 1, 37, 3, 28, 45, 10, 18, 24, 45, 17, 50, 0, 6, 14, 37, 50, 30, 22, 40, 4, 48, 37, 45, 10, 19, 38, 3, 15, 38, 9, 48, 50, 35, 27},
{16, 41, 6, 14, 3, 13, 21, 49, 34, 10, 30, 39, 47, 17, 38, 1, 0, 49, 21, 16, 16, 16, 21, 39, 13, 26, 28, 4, 14, 6, 21, 28, 1, 18, 46, 14, 15, 37, 0, 34, 24, 37, 28, 30, 24, 38, 34, 9, 6, 32},
{40, 23, 20, 35, 18, 20, 26, 40, 44, 32, 48, 44, 9, 38, 14, 37, 49, 0, 11, 17, 44, 30, 15, 44, 43, 2, 43, 41, 2, 9, 43, 4, 38, 27, 40, 33, 16, 18, 35, 6, 28, 48, 24, 27, 34, 47, 7, 23, 42, 41},
{39, 28, 12, 48, 49, 16, 47, 32, 19, 36, 2, 28, 5, 42, 19, 3, 21, 11, 0, 8, 26, 49, 33, 35, 19, 16, 29, 29, 11, 13, 12, 45, 21, 9, 44, 6, 49, 1, 44, 17, 39, 31, 17, 6, 10, 48, 20, 29, 49, 27},
{33, 39, 17, 47, 46, 5, 15, 35, 18, 23, 37, 22, 50, 35, 37, 28, 16, 17, 8, 0, 17, 18, 23, 35, 20, 48, 20, 26, 4, 37, 40, 7, 31, 25, 20, 49, 9, 49, 18, 18, 36, 29, 15, 38, 10, 30, 37, 4, 28, 23},
{35, 35, 21, 41, 42, 22, 46, 40, 18, 10, 2, 45, 23, 17, 13, 45, 16, 44, 26, 17, 0, 29, 12, 17, 22, 34, 11, 21, 25, 36, 42, 37, 13, 40, 8, 44, 43, 25, 30, 32, 41, 37, 15, 47, 32, 2, 27, 4, 36, 42},
{9, 35, 33, 11, 26, 14, 21, 38, 26, 37, 21, 8, 25, 5, 37, 10, 16, 30, 49, 18, 29, 0, 18, 40, 34, 1, 50, 45, 35, 5, 36, 42, 22, 46, 27, 50, 22, 46, 10, 11, 26, 17, 7, 35, 17, 33, 16, 48, 22, 37},
{0, 47, 45, 8, 36, 45, 37, 40, 20, 27, 24, 42, 46, 34, 5, 18, 21, 15, 33, 23, 12, 18, 0, 36, 30, 24, 24, 41, 44, 34, 28, 45, 3, 10, 47, 16, 4, 38, 44, 5, 38, 18, 46, 49, 1, 32, 28, 40, 20, 26},
{26, 28, 23, 14, 22, 29, 15, 47, 15, 47, 16, 2, 6, 43, 30, 24, 39, 44, 35, 35, 17, 40, 36, 0, 18, 38, 6, 26, 28, 13, 34, 6, 13, 17, 17, 15, 26, 47, 4, 2, 24, 4, 20, 45, 15, 7, 45, 4, 19, 4},
{3, 0, 7, 47, 46, 16, 0, 26, 16, 47, 1, 33, 33, 32, 21, 45, 13, 43, 19, 20, 22, 34, 30, 18, 0, 26, 26, 29, 19, 4, 24, 25, 38, 13, 3, 13, 9, 1, 36, 7, 0, 30, 30, 10, 37, 11, 47, 48, 40, 35},
{18, 35, 3, 17, 32, 28, 6, 24, 2, 17, 26, 5, 42, 46, 2, 17, 26, 2, 16, 48, 34, 1, 24, 38, 26, 0, 22, 35, 47, 44, 47, 1, 26, 5, 3, 49, 21, 35, 41, 32, 8, 27, 12, 15, 45, 20, 17, 7, 5, 15},
{17, 17, 43, 15, 20, 22, 5, 12, 6, 38, 47, 27, 7, 5, 32, 50, 28, 43, 29, 20, 11, 50, 24, 6, 26, 22, 0, 36, 13, 28, 24, 41, 18, 3, 26, 48, 3, 27, 26, 14, 6, 41, 10, 0, 11, 13, 42, 40, 44, 15},
{41, 49, 26, 2, 35, 30, 5, 39, 48, 14, 47, 11, 7, 37, 49, 0, 4, 41, 29, 26, 21, 45, 41, 26, 29, 35, 36, 0, 30, 25, 13, 11, 20, 48, 40, 30, 15, 10, 13, 36, 26, 44, 5, 47, 34, 27, 18, 27, 16, 41},
{4, 19, 36, 26, 18, 29, 23, 30, 7, 38, 5, 29, 48, 47, 43, 6, 14, 2, 11, 4, 25, 35, 44, 28, 19, 47, 13, 30, 0, 30, 35, 44, 1, 13, 15, 23, 12, 2, 29, 1, 30, 8, 13, 29, 16, 9, 9, 20, 26, 36},
{19, 34, 27, 27, 29, 15, 2, 40, 4, 14, 1, 15, 23, 39, 38, 14, 6, 9, 13, 37, 36, 5, 34, 13, 4, 44, 28, 25, 30, 0, 31, 31, 20, 28, 18, 6, 12, 39, 14, 41, 15, 3, 13, 24, 5, 36, 41, 3, 6, 8},
{49, 11, 0, 42, 8, 37, 49, 13, 20, 39, 0, 3, 38, 40, 23, 37, 21, 43, 12, 40, 42, 36, 28, 34, 24, 47, 24, 13, 35, 31, 0, 46, 7, 21, 7, 47, 39, 25, 28, 28, 7, 42, 43, 33, 33, 41, 25, 7, 46, 36},
{35, 15, 25, 9, 41, 42, 1, 7, 15, 4, 39, 5, 46, 10, 35, 50, 28, 4, 45, 7, 37, 42, 45, 6, 25, 1, 41, 11, 44, 31, 46, 0, 38, 26, 16, 13, 30, 11, 35, 2, 26, 5, 16, 33, 39, 15, 1, 13, 22, 50},
{19, 12, 4, 26, 17, 30, 45, 6, 16, 12, 21, 5, 3, 46, 8, 30, 1, 38, 21, 31, 13, 22, 3, 13, 38, 26, 18, 20, 1, 20, 7, 38, 0, 15, 43, 40, 8, 21, 9, 33, 16, 18, 34, 26, 31, 49, 40, 23, 25, 1},
{24, 27, 46, 46, 7, 29, 43, 16, 3, 17, 16, 49, 41, 46, 3, 22, 18, 27, 9, 25, 40, 46, 10, 17, 13, 5, 3, 48, 13, 28, 21, 26, 15, 0, 5, 41, 4, 34, 5, 46, 38, 29, 35, 23, 32, 30, 2, 26, 16, 17},
{24, 6, 10, 41, 1, 1, 18, 27, 46, 43, 13, 44, 15, 1, 30, 40, 46, 40, 44, 20, 8, 27, 47, 17, 3, 3, 26, 40, 15, 18, 7, 16, 43, 5, 0, 43, 5, 43, 40, 32, 37, 44, 25, 36, 29, 24, 28, 23, 7, 13},
{40, 23, 28, 46, 37, 33, 44, 13, 9, 18, 8, 47, 42, 45, 40, 4, 14, 33, 6, 49, 44, 50, 16, 15, 13, 49, 48, 30, 23, 6, 47, 13, 40, 41, 43, 0, 22, 43, 14, 16, 36, 39, 19, 5, 34, 43, 39, 43, 15, 23},
{25, 13, 15, 9, 13, 43, 24, 1, 22, 22, 22, 19, 28, 41, 2, 48, 15, 16, 49, 9, 43, 22, 4, 26, 9, 21, 3, 15, 12, 12, 39, 30, 8, 4, 5, 22, 0, 35, 34, 43, 2, 42, 17, 46, 38, 31, 15, 29, 26, 44},
{16, 11, 11, 10, 18, 3, 29, 50, 42, 3, 24, 36, 15, 34, 4, 37, 37, 18, 1, 49, 25, 46, 38, 47, 1, 35, 27, 10, 2, 39, 25, 11, 21, 34, 43, 43, 35, 0, 19, 35, 27, 9, 9, 24, 43, 15, 25, 8, 24, 14},
{19, 46, 11, 36, 26, 20, 40, 6, 49, 50, 50, 1, 18, 23, 14, 45, 0, 35, 44, 18, 30, 10, 44, 4, 36, 41, 26, 13, 29, 14, 28, 35, 9, 5, 40, 14, 34, 19, 0, 22, 24, 28, 7, 4, 17, 6, 26, 11, 11, 17},
{50, 6, 3, 1, 0, 2, 24, 7, 14, 18, 2, 16, 40, 50, 0, 10, 34, 6, 17, 18, 32, 11, 5, 2, 7, 32, 14, 36, 1, 41, 28, 2, 33, 46, 32, 16, 43, 35, 22, 0, 10, 35, 50, 17, 28, 36, 4, 30, 49, 29},
{5, 22, 16, 17, 39, 49, 34, 27, 46, 38, 35, 26, 45, 27, 1, 19, 24, 28, 39, 36, 41, 26, 38, 24, 0, 8, 6, 26, 30, 15, 7, 26, 16, 38, 37, 36, 2, 27, 24, 10, 0, 1, 38, 45, 24, 17, 45, 49, 5, 24},
{5, 13, 3, 10, 29, 39, 17, 22, 15, 8, 5, 9, 3, 27, 40, 38, 37, 48, 31, 29, 37, 17, 18, 4, 30, 27, 41, 44, 8, 3, 42, 5, 18, 29, 44, 39, 42, 9, 28, 35, 1, 0, 32, 39, 20, 25, 49, 29, 14, 13},
{5, 28, 36, 49, 15, 50, 12, 31, 36, 37, 21, 21, 1, 19, 31, 3, 28, 24, 17, 15, 15, 7, 46, 20, 30, 12, 10, 5, 13, 13, 43, 16, 34, 35, 25, 19, 17, 9, 7, 50, 38, 32, 0, 3, 42, 31, 36, 30, 4, 40},
{49, 31, 7, 46, 14, 47, 37, 50, 30, 1, 9, 44, 43, 46, 41, 15, 30, 27, 6, 38, 47, 35, 49, 45, 10, 15, 0, 47, 29, 24, 33, 33, 26, 23, 36, 5, 46, 24, 4, 17, 45, 39, 3, 0, 25, 29, 23, 16, 36, 50},
{7, 2, 24, 2, 17, 10, 38, 44, 44, 33, 19, 13, 36, 1, 22, 38, 24, 34, 10, 10, 32, 17, 1, 15, 37, 45, 11, 34, 16, 5, 33, 39, 31, 32, 29, 34, 38, 43, 17, 28, 24, 20, 42, 25, 0, 48, 50, 46, 26, 31},
{26, 41, 13, 7, 37, 45, 22, 32, 30, 22, 27, 29, 24, 30, 0, 9, 38, 47, 48, 30, 2, 33, 32, 7, 11, 20, 13, 27, 9, 36, 41, 15, 49, 30, 24, 43, 31, 15, 6, 36, 17, 25, 31, 29, 48, 0, 46, 41, 13, 46},
{43, 41, 8, 7, 27, 15, 7, 44, 23, 30, 1, 2, 17, 46, 9, 48, 34, 7, 20, 37, 27, 16, 28, 45, 47, 17, 42, 18, 9, 41, 25, 1, 40, 2, 28, 39, 15, 25, 26, 4, 45, 49, 36, 23, 50, 46, 0, 8, 10, 10},
{21, 20, 35, 26, 47, 25, 17, 12, 40, 9, 1, 41, 8, 14, 3, 50, 9, 23, 29, 4, 4, 48, 40, 4, 48, 7, 40, 27, 20, 3, 7, 13, 23, 26, 23, 43, 29, 8, 11, 30, 49, 29, 30, 16, 46, 41, 8, 0, 46, 25},
{9, 24, 47, 5, 32, 45, 49, 25, 7, 14, 45, 44, 21, 30, 0, 35, 6, 42, 49, 28, 36, 22, 20, 19, 40, 5, 44, 16, 26, 6, 46, 22, 25, 16, 7, 15, 26, 24, 11, 49, 5, 14, 4, 36, 26, 13, 10, 46, 0, 1},
{43, 16, 41, 31, 30, 8, 34, 23, 13, 32, 3, 33, 30, 14, 37, 27, 32, 41, 27, 23, 42, 37, 26, 4, 35, 15, 15, 41, 36, 8, 36, 50, 1, 17, 13, 23, 44, 14, 17, 29, 24, 13, 40, 50, 31, 46, 10, 25, 1, 0}};


  execute(mat1, mat2, sim, res);
  return 1;



  //std::vector<std::vector<double>> mat1(n), mat2(n);
  std::uniform_real_distribution<double> unif(-1000000,1000000);
  std::default_random_engine re;

  for (int i = 0; i < n; i++)
  {
    mat1[i].resize(n);
    mat2[i].resize(n);

    for (int j = 0; j < n; j++)
    {
      if (i < j)
      {
        mat1[i][j] = unif(re);
        mat2[i][j] = unif(re);
      }
      else if (i > j)
      {
        mat1[i][j] = mat1[j][i];
        mat2[i][j] = mat2[j][i];
      }
      else
      {
        mat1[i][j] = 0;
        mat2[i][j] = 0;
      }
    }
  }

  ttk::Timer timer;
  for (int i = 0; i < 20; i++)
    execute(mat1, mat2, sim, res);

  this->printMsg("CompleteAll", 1, timer.getElapsedTime());
  return 1;
}

int ttk::DistanceMatrixDistorsion::execute(const std::vector<std::vector<double>> &highDistMatrix, const std::vector<std::vector<double>> &lowDistMatrix, double &distorsionValue, std::vector<double> &distorsionVerticesValues) const
{
  ttk::Timer timer;
  auto n = highDistMatrix.size();

#ifndef TTK_ENABLE_KAMIKAZE
  // print horizontal separator
  this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator
                                             // print input parameters in table format
  this->printMsg({
      {"#Threads", std::to_string(this->threadNumber_)},
      {"#Vertices", std::to_string(n)},
      });
  this->printMsg(ttk::debug::Separator::L1);

  if (lowDistMatrix.size() != n)
  {
    this->printErr(" Sizes mismatch: the high distance matrix has " + std::to_string(n) + " rows and the low distance matrix has " + std::to_string(lowDistMatrix.size()) + " rows\n.");
    return 0;
  }
#endif

  /*
  this->printMsg("\n La high matrice :\n\n");
  for (size_t i = 0; i < n; i++)
  {
    for (size_t j = 0; j < n; j++)
      this->printMsg("\t"+std::to_string(highDistMatrix[i][j]));
    this->printMsg("\n");
  }
  this->printMsg("\n Et la low matrice :\n\n");
  for (size_t i = 0; i < n; i++)
  {
    for (size_t j = 0; j < n; j++)
      this->printMsg("\t"+std::to_string(lowDistMatrix[i][j]));
    this->printMsg("\n");
  }
  */

#ifndef TTK_ENABLE_KAMIKAZE
  for (size_t i = 0; i < n; i++)
  {
    if (highDistMatrix[i].size() != n)
    {
      this->printErr(" Sizes mismatch: high distance matrix is not a square matrix: it has " + std::to_string(n) + " rows and  row " + std::to_string(i) + " has " + std::to_string(highDistMatrix[i].size()) + " elements.\n");
      return 1;
    }
    if (lowDistMatrix[i].size() != n)
    {
      this->printErr(" Sizes mismatch: low distance matrix is not a square matrix: it has " + std::to_string(n) + " rows and  row " + std::to_string(i) + " has " + std::to_string(lowDistMatrix[i].size()) + "elements .\n");
      return 1;
    }
  }
#endif

  distorsionVerticesValues.resize(n);

  /* The computation, which is optimised for performance here, can be decomposed as follows:
   * compute for each (x,y) delta(x,y) = (dist_low(x,y)-dist_high(x,y))^2.
   * Compute maxi = maximum (delta(x,y)) over all (x,y).
   * Compute delta2(x,y) = 1-delta(x,y) for each (x,y).
   * The sim value is the mean of the n^2 values of delta2(x,y).
   * The distorsion for a vertex x0 is the mean of delta2(x0,y) over all y's.
  */

  double maxi = 0;
  /*std::vector<std::vector<double>> deltaBis(lowDistMatrix.size());

  for (size_t i = 0; i < n; i++)
    deltaBis[i].resize(n);*/

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_) reduction(max:maxi) schedule(dynamic)
#endif // TTK_ENABLE_OPENMP

  for (size_t i = 0; i < n; i++)
  {
    for (size_t j = i+1; j < n; j++) // plus lent si seulement calcul triangle supérieur
    {
      double diff = lowDistMatrix[i][j] - highDistMatrix[i][j];
      maxi = std::max(maxi, diff*diff);
      //deltaBis[i][j] = diff*diff;
      //maxi = std::max(maxi, deltaBis[i][j]);
    }
  }
  if (maxi < 1e-8) // TODO nbbits
    maxi = 1;
  //TODO et si maxi ~= 0 ?
  double totalSum = 0;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_), reduction(+:totalSum)
#endif // TTK_ENABLE_OPENMP

  //TODO attention parallèle sommes flottantes
  for (size_t i = 0; i < n; i++)
      //deltaBis[i][j] = diff*diff;
  {
    double sum = 0;
    for (size_t j = 0; j < n; j++)
    {
      double diff = lowDistMatrix[i][j] - highDistMatrix[i][j];
      double diff2 = diff*diff;
      //deltaBis[i][j] = 1-(diff2/maxi);
      //sum += deltaBis[i][j];
      sum += diff2;
    }
    double sumHarmonized = sum/maxi;
    distorsionVerticesValues[i] = 1-sumHarmonized/n;
    totalSum += 1-sumHarmonized/n;
  }

  distorsionValue = totalSum/n;
#ifndef TTK_ENABLE_KAMIKAZE
  this->printMsg("Size of output in ttk/base = " + std::to_string(distorsionVerticesValues.size()) + "\n");

  this->printMsg("Computed distorsion value: " + std::to_string(distorsionValue));
  this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator
  this->printMsg("Complete", 1, timer.getElapsedTime());
  this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
#endif
  return 0;
}
