#ifndef TJMODEL_SRC_MYUTIL_H
#define TJMODEL_SRC_MYUTIL_H

#include <stdlib.h>
#include <vector>

size_t GetNumofMps();

void Show(std::vector<size_t> v);

bool ParserBondDimension(int, char *[], std::vector<size_t> &);


//below is useless for t-J model
std::vector<size_t> GenerateDirectStateLabel(const size_t ly,
                                             const size_t lx,
                                             const size_t num_hole,
                                             const size_t ky_int);

std::vector<std::vector<size_t>> GenAllOrderedMomentumPairs(const size_t ly);

std::vector<std::vector<size_t>> GenAllMomentumPairs(const size_t ly);

std::vector<std::vector<size_t>> GenAllMomentumPairs3(const size_t ly);

bool ParserMeasureSite(const int argc, char *argv[],
                       size_t &start,
                       size_t &end);

std::vector<double> linspace(double start, double stop, int num_points);

#endif //TJMODEL_SRC_MYUTIL_H