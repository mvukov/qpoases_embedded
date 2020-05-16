# Copyright 2020 Milan Vukov. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import sys

import casadi
import numpy

import common

CHAIN_X0 = [0, 0, 0]


def get_chain_mass_ode(num_masses):
  x = casadi.SX.sym('x', 3, num_masses)
  dot_x = casadi.SX.sym('dot_x', 3, num_masses)
  x_end = casadi.SX.sym('x_end', 3)
  u = casadi.SX.sym('u', 3)

  y = casadi.veccat(x, dot_x, x_end)

  l = 0.033
  d = 1.0
  m = 0.03
  g = [0, 0, -9.81]

  def get_f(x_i, x_ip1):
    x_diff = x_ip1 - x_i
    return d * (1 - l / casadi.norm_2(x_diff)) * x_diff

  ddot_x = casadi.SX.zeros(3, num_masses)
  for i in range(0, num_masses):
    if i == 0:
      x_real_end = x[:, i + 1] if num_masses > 1 else x_end
      ddot_x[:, i] = 1. / m * (get_f(x[:, i], x_real_end) -
                               get_f(CHAIN_X0, x[:, i])) + g
    elif i < num_masses - 1:
      ddot_x[:, i] = 1. / m * (get_f(x[:, i], x[:, i + 1]) -
                               get_f(x[:, i - 1], x[:, i])) + g
    else:
      ddot_x[:, i] = 1. / m * (get_f(x[:, i], x_end) -
                               get_f(x[:, i - 1], x[:, i])) + g

  dot_y = casadi.veccat(dot_x, ddot_x, u)

  assert y.shape == dot_y.shape

  return u, y, dot_y


def get_equilibrium_states(num_masses, y, dot_y, x_end_value):
  assert y.shape == dot_y.shape
  assert y.shape[0] == 3 * (2 * num_masses + 1)

  x = y[0:3 * num_masses]
  x_end = y[-3:]
  eq_f = casadi.Function('eq_f', [x, x_end],
                         [dot_y[3 * num_masses:6 * num_masses]])
  solver = casadi.rootfinder('solver', 'newton', eq_f, {'abstol': 1e-10})

  x_end_value = numpy.asarray(x_end_value).flatten()

  assert x_end_value[2] == CHAIN_X0[2]
  x_initial = numpy.linspace(CHAIN_X0, x_end_value, 2 + num_masses).T
  x_initial = x_initial[:, 1:-1].T.flatten()
  eq_solution = solver(x_initial, x_end_value)

  eq_state = numpy.zeros(y.shape)
  eq_state[0:3 * num_masses] = eq_solution
  eq_state[-3:] = x_end_value.reshape((-1, 1))

  return eq_state


def solve_chain_mass_nmpc_qp(num_masses, num_intervals):
  time_step = 0.1

  u, y, dot_y = get_chain_mass_ode(num_masses)

  init_state = get_equilibrium_states(num_masses, y, dot_y, [0, 1, 0])
  ref_state = get_equilibrium_states(num_masses, y, dot_y, [1, 0, 0])

  # NMPC formulation:

  num_states = y.shape[0]
  num_controls = 3

  weights_states = numpy.ones(y.shape)
  weights_states[0:3 * num_masses] = 1e2
  weights_states[3 * num_masses:6 * num_masses] = 1e0
  weights_states[-3:] = 1e2

  weights_controls = numpy.ones((num_controls, 1)) * 1e-2

  # states 1..N
  states = casadi.SX.sym('states', num_states, num_intervals)
  # controls 0..N-1
  controls = casadi.SX.sym('controls', num_controls, num_intervals)

  weights_states_sqrt = numpy.sqrt(weights_states)
  weights_controls_sqrt = numpy.sqrt(weights_controls)

  objective_residuals = []
  ref_control = numpy.zeros(u.shape)
  for node in range(num_intervals):
    objective_residuals.append(
        (states[:, node] - ref_state) * weights_states_sqrt)
    objective_residuals.append(
        (controls[:, node] - ref_control) * weights_controls_sqrt)

  objective_residuals = casadi.veccat(*objective_residuals)

  ode = casadi.Function('ode', [u, y], [dot_y])

  rk4_k1 = ode(u, y)
  rk4_k2 = ode(u, y + time_step / 2.0 * rk4_k1)
  rk4_k3 = ode(u, y + time_step / 2.0 * rk4_k2)
  rk4_k4 = ode(u, y + time_step * rk4_k3)
  final_y = y + time_step / 6.0 * (rk4_k1 + 2 * rk4_k2 + 2 * rk4_k3 + rk4_k4)
  integrate = casadi.Function('integrate', [u, y], [final_y])

  states_for_integration = casadi.horzcat(init_state, states[:, :-1])
  integrated_states = integrate.map(num_intervals)(controls,
                                                   states_for_integration)
  equality_constraints = states - integrated_states
  equality_constraints = casadi.veccat(equality_constraints)

  # Prepare and condense the underlying QP.
  states_vec = casadi.veccat(states)
  controls_vec = casadi.veccat(controls)

  jac_obj_residuals_wrt_states = casadi.jacobian(objective_residuals,
                                                 states_vec)
  jac_obj_residuals_wrt_controls = casadi.jacobian(objective_residuals,
                                                   controls_vec)

  jac_eq_constraints_wrt_states = casadi.jacobian(equality_constraints,
                                                  states_vec)
  jac_eq_constraints_wrt_controls = casadi.jacobian(equality_constraints,
                                                    controls_vec)

  qp_h = jac_obj_residuals_wrt_controls.T @ jac_obj_residuals_wrt_controls

  delta_x_contrib = casadi.solve(jac_eq_constraints_wrt_states,
                                 jac_eq_constraints_wrt_controls)
  delta_x_qp_contrib = -jac_obj_residuals_wrt_states @ delta_x_contrib

  qp_h += delta_x_qp_contrib.T @ delta_x_qp_contrib

  qp_g = (jac_obj_residuals_wrt_controls +
          delta_x_qp_contrib).T @ objective_residuals

  states_lbx = numpy.zeros(y.shape)
  states_ubx = numpy.zeros(y.shape)

  states_lbx.fill(-numpy.inf)
  for m in range(num_masses):
    states_lbx[3 * m + 1] = -0.01
  states_ubx.fill(numpy.inf)

  qp_lb_a = []
  qp_ub_a = []

  for _ in range(num_intervals):
    qp_lb_a.append(states_lbx)
    qp_ub_a.append(states_ubx)
  qp_lb_a = numpy.concatenate(qp_lb_a, axis=0)
  qp_ub_a = numpy.concatenate(qp_ub_a, axis=0)

  init_states = numpy.concatenate([init_state for _ in range(num_intervals)],
                                  axis=0)
  delta_x_bound_contrib = casadi.solve(jac_eq_constraints_wrt_states,
                                       equality_constraints) + init_states

  qp_lb_a -= delta_x_bound_contrib
  qp_ub_a -= delta_x_bound_contrib

  qp_a = -delta_x_contrib

  qp_fcn = casadi.Function('qp_h_fcn', [controls_vec, states_vec],
                           [qp_h, qp_g, qp_a, qp_lb_a, qp_ub_a])

  init_controls = numpy.zeros(controls_vec.shape)
  qp_h_eval, qp_g_eval, qp_a_eval, qp_lb_a_eval, qp_ub_a_eval = qp_fcn(
      init_controls, init_states)

  qp_lbx = -numpy.ones(qp_g.shape)
  qp_ubx = numpy.ones(qp_g.shape)

  # Reduce the number of rows of the A-matrix to the minimum.
  qp_a_indices = []
  for el in range(qp_lb_a.shape[0]):
    if qp_lb_a_eval[el] <= -numpy.inf and qp_ub_a_eval[el] >= numpy.inf:
      continue
    qp_a_indices.append(el)

  qp_lb_a_eval = qp_lb_a_eval[qp_a_indices]
  qp_ub_a_eval = qp_ub_a_eval[qp_a_indices]
  qp_a_eval = qp_a_eval[qp_a_indices, :]

  qp_x = casadi.SX.sym('qp_x', *qp_g.shape)
  qp = {
      'x': qp_x,
      'f': 0.5 * qp_x.T @ qp_h_eval @ qp_x + qp_x.T @ qp_g_eval,
      'g': casadi.densify(qp_a_eval @ qp_x)
  }

  qp_solver = casadi.qpsol('qp_solver', 'qpoases', qp)

  sol = qp_solver(lbx=qp_lbx, ubx=qp_ubx, lbg=qp_lb_a_eval, ubg=qp_ub_a_eval)

  x_opt = numpy.asarray(sol['x']).flatten()
  y_opt = numpy.asarray(-casadi.vertcat(sol['lam_x'], sol['lam_g']))
  f_opt = sol['f']

  num_variables = qp_x.shape[0]
  num_constraints = qp_a_eval.shape[0]

  qp_h_flat = numpy.asarray(qp_h_eval).flatten()
  qp_g_flat = numpy.asarray(qp_g_eval).flatten()
  qp_a_flat = numpy.asarray(qp_a_eval).flatten()
  qp_lb_a_flat = numpy.asarray(qp_lb_a_eval).flatten()
  qp_ub_a_flat = numpy.asarray(qp_ub_a_eval).flatten()

  return (num_variables, num_constraints, qp_h_flat, qp_g_flat, qp_a_flat,
          qp_lbx, qp_ubx, qp_lb_a_flat, qp_ub_a_flat, x_opt, y_opt, f_opt)


def main():
  test_data_vector = [
      solve_chain_mass_nmpc_qp(num_masses=num_masses,
                               num_intervals=num_intervals)
      for num_masses in range(1, 4)
      for num_intervals in [10, 20]
  ]

  header = common.get_test_data_header(test_data_vector)
  common.write_to_file(header, sys.argv[1])


if __name__ == '__main__':
  main()
