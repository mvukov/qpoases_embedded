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
import os

import numpy


def get_test_data_string(test_data):
  numpy.set_printoptions(precision=16)
  test_data_str = []
  for data in test_data:
    if isinstance(data, numpy.ndarray):
      flat_data = data.flatten()
      test_data_str.append('{{{}}}'.format(', '.join(
          [str(value) for value in flat_data])))
    else:
      test_data_str.append(str(data))
  return '{{{}}}'.format(', '.join(test_data_str))


def get_test_data_header(test_data_vector):
  test_data_vector_strings = ',\n'.join(
      [get_test_data_string(test_data) for test_data in test_data_vector])

  header = """
#include "qpoases_embedded/Constants.hpp"

namespace qpoases_embedded {{

struct QpTestData {{
  int num_variables;
  int num_constraints;
  std::vector<real_t> h;
  std::vector<real_t> g;
  std::vector<real_t> a;
  std::vector<real_t> lb;
  std::vector<real_t> ub;
  std::vector<real_t> lba;
  std::vector<real_t> uba;
  std::vector<real_t> x_opt;
  std::vector<real_t> y_opt;
  real_t f_opt;
}};

static constexpr real_t inf = INFTY;

const std::vector<QpTestData> qp_test_data_vectors = {{
{test_data_vector_strings}
}};

}} // namespace qpoases_embedded
""".format(test_data_vector_strings=test_data_vector_strings)
  return header


def write_to_file(data, filename):
  file_dir = os.path.dirname(filename)
  if file_dir and not os.path.exists(file_dir):
    os.makedirs(file_dir)
  with open(filename, 'w') as stream:
    stream.write(data)
