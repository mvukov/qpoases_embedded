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
#include <algorithm>
#include <optional>
#include <vector>

#include <benchmark/benchmark.h>

#include <qpoases_embedded/MessageHandling.hpp>
#include <qpoases_embedded/QProblem.hpp>
#include <qpoases_embedded/Utils.hpp>

#include "./qp_test_data.h"

#include "qpoases_embedded/chain_mass_nmpc_benchmark_data.h"

namespace qpoases_embedded {

void PrepareBenchmarkArguments(benchmark::internal::Benchmark* b) {
  for (const auto& test_data : qp_test_data_vectors) {
    b->Args({test_data.num_variables, test_data.num_constraints});
  }
}

std::optional<QProblem> GetQpSolver(benchmark::State& state) {
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

void BenchmarkChainMassNmpcQp(benchmark::State& state) {
  auto qp_solver = GetQpSolver(state);
  if (!qp_solver) {
    state.SkipWithError("Failed to find QP test data!");
    return;
  }

  static constexpr size_t kNumSolvers = 1000;
  std::vector<QProblem> qp_solvers(kNumSolvers, *qp_solver);

  size_t index = 0;
  for (auto _ : state) {
    state.PauseTiming();
    auto& current_qp_solver = qp_solvers.at(index);
    std::vector<double>* qp_x = current_qp_solver.getMutableX();
    std::vector<double>* qp_y = current_qp_solver.getMutableY();
    std::fill(qp_x->begin(), qp_x->end(), 0);
    std::fill(qp_y->begin(), qp_y->end(), 0);
    state.ResumeTiming();

    int nWSR = 100;
    const auto success = current_qp_solver.init(nWSR, {});
    if (success != SUCCESSFUL_RETURN) {
      state.SkipWithError("Failed to solve QP!");
      break;
    }
    index = (index + 1) % kNumSolvers;
  }
}

BENCHMARK(BenchmarkChainMassNmpcQp)
    ->Apply(PrepareBenchmarkArguments)
    ->Unit(::benchmark::kMillisecond)
    ->MinTime(2.0);

void BenchmarkChainMassNmpcQpCached(benchmark::State& state) {
  auto qp_solver = GetQpSolver(state);
  if (!qp_solver) {
    state.SkipWithError("Failed to find QP test data!");
    return;
  }

  std::vector<double>* qp_x = qp_solver->getMutableX();
  std::vector<double>* qp_y = qp_solver->getMutableY();

  for (auto _ : state) {
    state.PauseTiming();
    std::fill(qp_x->begin(), qp_x->end(), 0);
    std::fill(qp_y->begin(), qp_y->end(), 0);
    state.ResumeTiming();

    int nWSR = 100;
    const auto success = qp_solver->init(nWSR, {});
    if (success != SUCCESSFUL_RETURN) {
      state.SkipWithError("Failed to solve QP!");
      break;
    }
  }
}

BENCHMARK(BenchmarkChainMassNmpcQpCached)
    ->Apply(PrepareBenchmarkArguments)
    ->Unit(::benchmark::kMillisecond)
    ->MinTime(2.0);

}  // namespace qpoases_embedded

BENCHMARK_MAIN();
