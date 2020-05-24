/*
 *  This file is part of qpOASES.
 *
 *  qpOASES -- An Implementation of the Online Active Set Strategy.
 *  Copyright (C) 2020 Milan Vukov.
 *
 *  qpOASES is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  qpOASES is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with qpOASES; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301
 *  USA
 *
 */
#ifndef BENCHMARK_HELPERS_H_
#define BENCHMARK_HELPERS_H_

#include <algorithm>
#include <optional>

#include <benchmark/benchmark.h>

namespace qpoases_embedded {

// qp_test_data_vectors: vector of test data that is automatically generated.
// This vector is provided in an auto-generated file.

void PrepareBenchmarkArguments(benchmark::internal::Benchmark* b) {
  for (const auto& test_data : qp_test_data_vectors) {
    b->Args({test_data.num_variables, test_data.num_constraints});
  }
}

std::optional<QProblem> GetQpSolver(const benchmark::State& state) {
  const auto test_data =
      std::find_if(qp_test_data_vectors.begin(), qp_test_data_vectors.end(),
                   [&state](const auto& t) {
                     return t.num_variables == state.range(0) &&
                            t.num_constraints == state.range(1);
                   });
  if (test_data == qp_test_data_vectors.end()) {
    return {};
  }

  QProblem qp_solver(test_data->num_variables, test_data->num_constraints);
  *qp_solver.getMutableH() = test_data->h;
  *qp_solver.getMutableG() = test_data->g;
  *qp_solver.getMutableA() = test_data->a;
  *qp_solver.getMutableLb() = test_data->lb;
  *qp_solver.getMutableUb() = test_data->ub;
  *qp_solver.getMutableLbA() = test_data->lba;
  *qp_solver.getMutableUbA() = test_data->uba;

  return qp_solver;
}

}  // namespace qpoases_embedded

#endif  // BENCHMARK_HELPERS_H_
