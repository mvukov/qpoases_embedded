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
#include <vector>

#include <qpoases_embedded/MessageHandling.hpp>
#include <qpoases_embedded/QProblem.hpp>

#include "qpoases_embedded/hanging_chain_test_data.h"

#include "./benchmark_helpers.h"

namespace qpoases_embedded {

void BenchmarkHangingChainQpCachedIteration0(benchmark::State& state) {
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
    const auto success =
        qp_solver->init(nWSR, [](int, real_t, int, int, int, SubjectToStatus,
                                 bool) { return false; });
    if (success != RET_USER_ABORT_REQUESTED && success != SUCCESSFUL_RETURN) {
      state.SkipWithError("Failed to solve QP!");
      break;
    }
  }
}

BENCHMARK(BenchmarkHangingChainQpCachedIteration0)
    ->Apply(PrepareBenchmarkArguments)
    ->Unit(::benchmark::kMillisecond)
    ->MinTime(1.0);

void BenchmarkHangingChainQpCachedIteration1(benchmark::State& state) {
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

    int nWSR = 100;
    const auto success = qp_solver->init(
        nWSR,
        [&state](int iteration, real_t, int, int, int, SubjectToStatus, bool) {
          if (iteration == 0) {
            state.ResumeTiming();
          } else if (iteration == 1) {
            return false;
          }
          return true;
        });
    if (success != RET_USER_ABORT_REQUESTED && success != SUCCESSFUL_RETURN) {
      state.SkipWithError("Failed to solve QP!");
      break;
    }
  }
}

BENCHMARK(BenchmarkHangingChainQpCachedIteration1)
    ->Apply(PrepareBenchmarkArguments)
    ->Unit(::benchmark::kMillisecond)
    ->MinTime(1.0);

}  // namespace qpoases_embedded

BENCHMARK_MAIN();
