// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef NETKET_PYGRAPH_HPP
#define NETKET_PYGRAPH_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include "abstract_graph.hpp"

namespace py = pybind11;

namespace netket {

namespace detail {
/// Given a Python iterable, returns its length if it is known or 0 otherwise.
/// This can safely be used to preallocate storage on the C++ side as calls to
/// `std::vector<T>::reserve(0)` are basically noops.
inline std::size_t LengthHint(py::iterable xs) {
  auto iterator = xs.attr("__iter__")();
  if (py::hasattr(iterator, "__length_hint__")) {
    auto const n = iterator.attr("__length_hint__")().cast<long>();
    assert(n >= 0 && "Bug in Python/pybind11??");
    return static_cast<std::size_t>(n);
  }
  return 0;
}
inline std::size_t LengthHint(py::iterator x) {
  if (py::hasattr(x, "__length_hint__")) {
    auto const n = x.attr("__length_hint__")().cast<long>();
    assert(n >= 0 && "Bug in Python/pybind11??");
    return static_cast<std::size_t>(n);
  }
  return 0;
}

/// Correctly orders site indices and constructs an edge.
// TODO(twesterhout): Should we throw when `x == y`? I.e. edge from a node to
// itself is a questionable concept.
inline AbstractGraph::Edge MakeEdge(int const x, int const y) noexcept {
  using Edge = AbstractGraph::Edge;
  return (x < y) ? Edge{{x, y}} : Edge{{y, x}};
}

/// Converts a Python iterable to a list of edges. An exception is thrown if the
/// input iterable contains duplicate edges.
///
/// \postcondition For each edge (i, j) we have i <= j.
/// \postcondition The returned list contains no duplicates.
inline std::vector<AbstractGraph::Edge> Iterable2Edges(py::iterator x) {
  using std::begin;
  using std::end;
  std::vector<AbstractGraph::Edge> edges;
  edges.reserve(LengthHint(x));

  while (x != py::iterator::sentinel()) {
    int i, j;
    std::tie(i, j) = x->template cast<std::tuple<int, int>>();
    edges.push_back(MakeEdge(i, j));
    ++x;
  }

  // NOTE(twesterhout): Yes, I know that this screws up the algorithmic
  // complexity, but it's fast enough to be unnoticeable for all practical
  // purposes.
  std::sort(begin(edges), end(edges));
  if (std::unique(begin(edges), end(edges)) != end(edges)) {
    throw InvalidInputError{"Edge list contains duplicates."};
  }
  return edges;
}

/// Converts a Python iterable to a `ColorMap`. An exception is thrown if the
/// input iterable contains duplicate edges.
///
/// \postcondition For each edge (i, j) we have i <= j.
inline AbstractGraph::ColorMap Iterable2ColorMap(py::iterator x) {
  AbstractGraph::ColorMap colors;
  colors.reserve(LengthHint(x));

  while (x != py::iterator::sentinel()) {
    int i, j, color;
    std::tie(i, j, color) = x->template cast<std::tuple<int, int, int>>();
    if (!colors.emplace(MakeEdge(i, j), color).second) {
      // Failed to insert an edge because it already exists
      throw InvalidInputError{"Edge list contains duplicates."};
    }
    ++x;
  }
  return colors;
}
}  // namespace detail
}  // namespace netket

#include "py_custom_graph.hpp"
#include "py_lattice.hpp"

namespace netket {

void AddGraphModule(py::module& m) {
  auto subm = m.def_submodule("graph");

  py::class_<AbstractGraph>(subm, "Graph")
      .def_property_readonly("n_sites", &AbstractGraph::Nsites,
                             R"EOF(
      int: The number of vertices in the graph.)EOF")
      .def_property_readonly(
          "edges",
          [](AbstractGraph const& x) {
            using vector_type =
                std::remove_reference<decltype(x.Edges())>::type;
            return vector_type{x.Edges()};
          },
          R"EOF(
      list: The graph edges.)EOF")
      .def_property_readonly("adjacency_list", &AbstractGraph::AdjacencyList,
                             R"EOF(
      list: The adjacency list of the graph where each node is
          represented by an integer in `[0, n_sites)`.)EOF")
      .def_property_readonly("is_bipartite", &AbstractGraph::IsBipartite,
                             R"EOF(
      bool: Whether the graph is bipartite.)EOF")
      .def_property_readonly("is_connected", &AbstractGraph::IsConnected,
                             R"EOF(
      bool: Whether the graph is connected.)EOF")
      .def_property_readonly("distances", &AbstractGraph::AllDistances,
                             R"EOF(
      list[list]: The distances between the nodes. The fact that some node
          may not be reachable from another is represented by -1.)EOF")
      .def_property_readonly("automorphisms", &AbstractGraph::SymmetryTable,
                             R"EOF(
      list[list]: The automorphisms of the graph,
          including translation symmetries only.)EOF")
      .def_static("Hypercube",
                  [](int length, int n_dim, bool pbc) -> Lattice {
                    std::vector<int> extent(n_dim, length);
                    std::vector<bool> pbc_vec(n_dim, pbc);
                    std::vector<std::vector<double>> basis_vectors(
                        n_dim, std::vector<double>(n_dim));
                    std::vector<std::vector<double>> atoms_coord(
                        1, std::vector<double>(n_dim));
                    for (int i = 0; i < n_dim; i++) {
                      atoms_coord[0][i] = 0;
                      for (int j = 0; j < n_dim; j++) {
                        if (i == j) {
                          basis_vectors[i][j] = 1;
                        } else {
                          basis_vectors[i][j] = 0;
                        }
                      }
                    }
                    return Lattice(basis_vectors, extent, pbc_vec, atoms_coord);
                  },
                  py::arg("length"), py::arg("n_dim"), py::arg("pbc") = true,
                  R"EOF(
            Member function constructing a hypercubic lattice of arbitrary dimension.

            Args:
                length: Side length of the hypercube.
                n_dim: Dimension of the hypercube.
                pbc: If ``True`` then the constructed hypercube
                    will have periodic boundary conditions, otherwise
                    open boundary conditions are imposed.
            )EOF");

  AddCustomGraph(subm);
  AddLattice(subm);

}  // namespace netket

}  // namespace netket

#endif
