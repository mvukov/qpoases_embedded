#!/usr/bin/env python3
#
# CasADi -- A symbolic framework for dynamic optimization.
# Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
#                         K.U. Leuven. All rights reserved.
# Copyright (C) 2011-2014 Greg Horn
#
# CasADi is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# CasADi is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with CasADi; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# This file is part of qpOASES.
#
# qpOASES -- An Implementation of the Online Active Set Strategy.
# Copyright (C) 2020 Milan Vukov.
#
# qpOASES is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# qpOASES is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with qpOASES; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301
# USA
"""
We want to model a chain attached to two supports and hanging in between.
Let us discretise it with N mass points connected by N - 1 springs. Each mass i
has position (y_i, z_i), i = 1,...,N. The equilibrium point of the system
minimises the potential energy.

The potential energy of each spring is
Vi = D_i / 2 * ((y_i - y_{i + 1})^2 + (z_i - z_{i + 1})^2)

The gravitational potential energy of each mass is
Vg_i = m_i * g0 * z_i

The total potential energy is thus given by:

Vchain(y, z) = 1/2 * sum{i = 1,..., N - 1} * D_i *
  ((y_i - y_{i + 1})^2 + (z_i - z_{i + 1})^2) + g0 * sum{i = 1,..., N} m_i * z_i

where y = [y_1,..., y_N] and z = [z_1,..., z_N]

We wish to minimize the following objective subject to ground constraints:

minimize {y, z} Vchain(y, z)
subject to: z_i >= zin
            z_i - 0.1*y_i >= 0.5
"""
import os
import sys

import casadi
import numpy


def solve_hanging_chain_qp(num_masses, use_contraints):
  m_i = 40.0 / num_masses
  D_i = 70.0 * num_masses
  g0 = 9.81
  zmin = 0.5  # ground

  x = []
  f = 0
  g = []

  lbx = []
  ubx = []
  lbg = []
  ubg = []

  y_start, z_start = -2, 1
  y_end, z_end = 2, 1

  # Loop over all chain elements
  for i in range(0, num_masses):
    # Previous point
    if i > 0:
      y_prev = y_i
      z_prev = z_i

    # Create variables for the (y_i, z_i) coordinates
    y_i = casadi.SX.sym('y_' + str(i))
    z_i = casadi.SX.sym('z_' + str(i))

    # Add to the list of variables
    x += [y_i, z_i]

    lbx += [-numpy.inf, zmin]
    ubx += [numpy.inf, numpy.inf]

    # Spring potential
    if i == 0:
      f += D_i / 2 * ((y_start - y_i)**2 + (z_start - z_i)**2)
    else:
      f += D_i / 2 * ((y_prev - y_i)**2 + (z_prev - z_i)**2)

    # Graviational potential
    f += g0 * m_i * z_i

    # Slanted ground constraints
    if use_contraints:
      g.append(z_i - 0.1 * y_i)
      lbg.append(0.5)
      ubg.append(numpy.inf)

  f += D_i / 2 * ((y_i - y_end)**2 + (z_i - z_end)**2)

  num_variables = len(x)
  num_constraints = len(g)

  x = casadi.vertcat(*x)
  g = casadi.vertcat(*g)
  qp = {'x': x, 'f': f, 'g': g}
  solver = casadi.qpsol('solver', 'qpoases', qp)

  gradient_f = casadi.gradient(f, x)
  h = casadi.jacobian(gradient_f, x, {'symmetric': True})
  c = casadi.substitute(gradient_f, x, casadi.SX.zeros(x.sparsity()))
  a = casadi.jacobian(g, x)

  prob = casadi.Function("qp_eval", [x], [h, c, a])
  qp_items = prob(casadi.SX.sym('eval_args', x.shape))
  h_flat = numpy.asarray(casadi.DM(qp_items[0]))
  c_flat = numpy.asarray(casadi.DM(qp_items[1]))
  a_flat = numpy.asarray(casadi.DM(qp_items[2]))

  sol = solver(lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg)

  x_opt = numpy.asarray(sol['x']).flatten()
  y_opt = numpy.asarray(-casadi.vertcat(sol['lam_x'], sol['lam_g']))
  f_opt = sol['f']

  return (num_variables, num_constraints, h_flat, c_flat, a_flat,
          numpy.asarray(lbx), numpy.asarray(ubx), numpy.asarray(lbg),
          numpy.asarray(ubg), x_opt, y_opt, f_opt)


def get_test_data_string(test_data):
  test_data_str = []
  for data in test_data:
    if isinstance(data, numpy.ndarray):
      flat_data = data.flatten()
      test_data_str.append('{{{}}}'.format(', '.join(
          [str(value) for value in flat_data])))
    else:
      test_data_str.append(str(data))
  return '{{{}}}'.format(', '.join(test_data_str))


def main():
  numpy.set_printoptions(precision=16)

  test_data_vector = [
      solve_hanging_chain_qp(num_masses=num_masses,
                             use_contraints=use_contraints)
      for num_masses in range(5, 100, 5)
      for use_contraints in [True, False]
  ]
  test_data_vector_strings = ',\n'.join(
      [get_test_data_string(test_data) for test_data in test_data_vector])

  header = """
const std::vector<QpTestData> qp_test_data_vectors = {{
{test_data_vector_strings}
}};
""".format(test_data_vector_strings=test_data_vector_strings)

  print(f'Header path: {sys.argv[1]}')

  header_path = sys.argv[1]
  header_dir = os.path.dirname(header_path)
  if header_dir and not os.path.exists(header_dir):
    os.makedirs(header_dir)
  with open(header_path, 'w') as stream:
    stream.write(header)


if __name__ == '__main__':
  main()
