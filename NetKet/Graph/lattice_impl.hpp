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

#ifndef NETKET_LATTICE_IMPL_HPP
#define NETKET_LATTICE_IMPL_HPP

#include <cmath>
#include <iostream>
#include <vector>
#include "Utils/array_search.hpp"
#include "Utils/messages.hpp"
#include "Utils/next_variation.hpp"

namespace netket {
// Constructor
Lattice::Lattice(std::vector<std::vector<double>> basis_vector,
                 std::vector<int> extent, std::vector<bool> pbc,
                 std::vector<std::vector<double>> V_atoms)
    : extent_(std::move(extent)),
      pbc_(std::move(pbc)),
      basis_vectors_(std::move(basis_vector)),
      atoms_coord_(std::move(V_atoms)) {
  ndim_ = basis_vectors_.size();

  // check if dimensions are consistents
  for (int i = 0; i < ndim_; i++) {
    if (basis_vectors_[i].size() != ndim_) {
      ErrorMessage() << "Each element of basis_vectors must have ndim"
                        "components. Dimensions are inconsistent. Exit.\n";
      exit(EXIT_FAILURE);
    }
  }

  if (extent_.size() != ndim_) {
    ErrorMessage() << "Extent must have ndim components. Dimension is "
                      "inconsistent. Exit.\n";
    exit(EXIT_FAILURE);
  }

  // default case is a single atom in the unit cell, located at the origin
  if (atoms_coord_.size() == 0) {
    atoms_coord_.resize(1, std::vector<double>(ndim_, 0));
  }

  natoms_ = atoms_coord_.size();

  for (int i = 0; i < natoms_; i++) {
    if (atoms_coord_[i].size() != ndim_) {
      ErrorMessage() << "Each element of atoms_coord must have ndim "
                        "components. Dimension is inconsistent. Exit.\n";
      exit(EXIT_FAILURE);
    }
  }

  // default case for pbc is true
  if (pbc_.size() == 0) {
    pbc_.resize(ndim_, true);
  }

  if (pbc_.size() != ndim_) {
    ErrorMessage() << "Pbc must have ndim components. Dimension is "
                      "inconsistent. Exit.\n";
    exit(EXIT_FAILURE);
  }

  // PATHOLOGIC 1: pbc(i) for extent(i)<=2 + extent(i)=0)
  for (int i = 0; i < ndim_; i++) {
    if (extent_[i] <= 2 and pbc_[i] == true) {
      ErrorMessage() << "PBC are not allowed when extent<2. Exit.\n";
      exit(EXIT_FAILURE);
    }
    if (extent_[i] == 0) {
      ErrorMessage() << "Extent components must be >0. Exit.\n";
      exit(EXIT_FAILURE);
    }
  }

  nlatticesites_ = 1;
  for (int k = 0; k < ndim_; k++) nlatticesites_ *= extent_[k];

  // PATHOLOGIC 2: nlatticesites_=0,1
  if (nlatticesites_ <= 1) {
    ErrorMessage()
        << "A well-defined lattice must have at least 2 sites. Exit.\n";
    exit(EXIT_FAILURE);
  }

  nsites_ = nlatticesites_ * natoms_;

  for (int k = 0; k < natoms_; k++) {
    for (int i = 0; i < nlatticesites_; i++) {
      std::vector<int> n = Site2Vector(i);  // from site i to coefficents
      std::vector<double> R =
          Vector2Coord(n, k);  // initialize the R vector (coord of one site)
      R_.push_back(R);
    }
  }

  edges_ = BuildEdges();

  symmetrytable_ = BuildSymmetryTable();

  colors_.reserve(nlatticesites_);
  for (auto const &edge : edges_) {
    auto success = colors_.emplace(edge, 0).second;
    static_cast<void>(success);  // Make everyone happy in the NDEBUG case
    assert(success && "There should be no duplicate edges");
  }

  is_connected_ = ComputeConnected();
  is_bipartite_ = ComputeBipartite();

};  // namespace netket
Lattice::~Lattice(){};

// Get private member
int Lattice::Ndim() const noexcept { return ndim_; };
int Lattice::Nsites() const noexcept { return nsites_; };
int Lattice::Size() const noexcept { return Nsites(); };

std::vector<std::vector<double>> Lattice::Coordinates() const noexcept {
  std::vector<std::vector<double>> result;
  for (int k = 0; k < nsites_; k++) result.push_back(R_[k]);
  return result;
}
std::vector<Lattice::Edge> const &Lattice::Edges() const noexcept {
  return edges_;
};
std::vector<std::vector<int>> Lattice::SymmetryTable() const noexcept {
  return symmetrytable_;
};

// Returns map of the edge and its respective color
const AbstractGraph::ColorMap &Lattice::EdgeColors() const noexcept {
  return colors_;
}

bool Lattice::ComputeConnected() const {
  const int start = 0;  // arbitrary node
  int nvisited = 0;
  BreadthFirstSearch(start, [&nvisited](int, int) { ++nvisited; });
  return nvisited == Nsites();
}

bool Lattice::ComputeBipartite() const {
  bool is_bipartite = true;
  const int start = 0;  // arbitrary node
  std::vector<int> colors(Nsites(), -1);
  const auto adjacency_list =
      AdjacencyList();  // implicit expression can't have
                        // access to the Lattice function
  if (is_connected_) {
    BreadthFirstSearch(
        start, [&colors, &adjacency_list, &is_bipartite](int node, int) {
          if (node == start) colors[node] = 1;
          for (std::size_t j = 0; j < adjacency_list[node].size(); j++) {
            if (!is_bipartite) break;
            if (colors[adjacency_list[node][j]] == -1) {
              colors[adjacency_list[node][j]] = 1 - colors[node];
            } else if (colors[adjacency_list[node][j]] == colors[node]) {
              is_bipartite = false;
            }
          }
        });
  } else {
    BreadthFirstSearch(
        [&colors, &adjacency_list, &is_bipartite](int node, int, int) {
          if (node == start) colors[node] = 1;
          for (std::size_t j = 0; j < adjacency_list[node].size(); j++) {
            if (!is_bipartite) break;
            if (colors[adjacency_list[node][j]] == -1) {
              colors[adjacency_list[node][j]] = 1 - colors[node];
            } else if (colors[adjacency_list[node][j]] == colors[node]) {
              is_bipartite = false;
            }
          }
        });
  }
  return is_bipartite;
}

// Lattice representations (site = k, vector = n_i, coord = cartesian
// coordinates)
std::vector<int> Lattice::Site2Vector(int i) const {
  int ndim = extent_.size();
  std::vector<int> result(ndim, 0);
  int ip;
  if (i < nlatticesites_)
    ip = i;
  else
    ip = i % nlatticesites_;
  int k = ndim - 1;
  while (ip > 0) {
    result[k] = ip % extent_[k];
    ip /= extent_[k];
    k--;
  }
  return result;
};
std::vector<double> Lattice::Vector2Coord(const std::vector<int> &n,
                                          int iatom) const {
  std::vector<double> R(ndim_, 0);
  for (int j = 0; j < ndim_; j++) {
    R[j] = atoms_coord_[iatom][j];
    for (int k = 0; k < ndim_; k++) R[j] += n[k] * basis_vectors_[k][j];
  }

  return R;
};
std::vector<double> Lattice::Site2Coord(int k) const { return R_[k]; }
int Lattice::Vector2Site(const std::vector<int> &n) const {
  int k = 0;

  for (int i = 0; i < Ndim(); i++) {
    int base = 1;
    for (int j = i + 1; j < Ndim(); j++) base *= extent_[j];
    k += n[i] * base;
  }

  return k;
}

// Neighbours
std::vector<std::vector<int>> Lattice::PossibleLatticeNeighbours() const {
  const int nsymbols = 3;
  const int max = nsymbols - 1;
  std::vector<int> b(ndim_, 0);
  std::vector<int> result_vector;
  std::vector<std::vector<int>> result_matrix;
  std::vector<int> row(ndim_, 0);
  do {
    for (std::vector<int>::const_iterator ib = b.begin(); ib != b.end(); ++ib)
      result_vector.push_back(*ib);
  } while (next_variation(b.begin(), b.end(), max));

  for (std::size_t i = 0; i < result_vector.size() / ndim_; i++) {
    for (int j = 0; j < ndim_; j++) {
      row[j] = result_vector[i * ndim_ + j] -
               1;  //-1 because I want (-1,0,1) as symbols instead of (0,1,2)
    }
    result_matrix.push_back(row);
  }
  return result_matrix;
};
std::vector<double> Lattice::NeighboursSquaredDistance(
    const std::vector<std::vector<int>> &neighbours_matrix, int iatom) const {
  std::vector<double> distance;
  for (std::size_t row = 0; row < neighbours_matrix.size(); row++) {
    std::vector<int> n(ndim_, 0);
    for (int j = 0; j < ndim_; j++) n[j] = neighbours_matrix[row][j];
    for (int jatom = 0; jatom < natoms_; jatom++)
      distance.push_back(
          GetSquaredDistance(Vector2Coord(n, jatom), atoms_coord_[iatom]));
  }

  return distance;
};
std::vector<std::vector<int>> Lattice::LatticeNeighbours(int iatom) const {
  std::vector<std::vector<int>> neighbours_matrix_in;
  std::vector<std::vector<int>> result;
  std::vector<double> distance;
  double min_distance;

  neighbours_matrix_in = PossibleLatticeNeighbours();

  distance = NeighboursSquaredDistance(neighbours_matrix_in, iatom);

  min_distance = *min_nonzero_elem(distance.begin(), distance.end());

  for (std::size_t i = 0; i < distance.size(); i++) {
    if (RelativelyEqual(distance[i], min_distance, 1.0e-6)) {
      std::vector<int> temp;
      temp = neighbours_matrix_in[i / natoms_];
      temp.push_back(i % natoms_);
      result.push_back(temp);
    }
  }

  return result;
}

std::vector<int> Lattice::FindNeighbours(int k, int iatom) const {
  std::vector<int> result;
  std::vector<int> n = Site2Vector(k);
  std::vector<int> single_neighbour;

  std::vector<std::vector<int>> lattice_neighbours_vector =
      LatticeNeighbours(iatom);

  for (std::size_t i = 0; i < lattice_neighbours_vector.size(); i++) {
    std::vector<int> single_neighbour_vector;

    bool flag = true;
    for (int j = 0; j < Ndim(); j++) {
      int new_index = n[j] + lattice_neighbours_vector[i][j];
      if (pbc_[j]) {
        if (new_index < 0)
          new_index = extent_[j] - 1;
        else if (new_index >= extent_[j])
          new_index = 0;
        single_neighbour_vector.push_back(new_index);

      } else if (new_index >= 0 and new_index < extent_[j])
        single_neighbour_vector.push_back(new_index);

      else
        flag = false;  // if new_index refers to a border site, don't add this
                       // site to the neighbours
    }
    if (flag) {
      single_neighbour.push_back(Vector2Site(single_neighbour_vector) +
                                 lattice_neighbours_vector[i][Ndim()] *
                                     nlatticesites_);
    }
  }

  return single_neighbour;
}

// Edges and AdjacencyList
std::vector<Lattice::Edge> Lattice::BuildEdges() const {
  std::vector<Lattice::Edge> edge_vector;
  for (int k = 0; k < nlatticesites_; k++) {
    for (int iatom = 0; iatom < natoms_; iatom++) {
      std::vector<int> neighbours;

      neighbours = Lattice::FindNeighbours(k, iatom);

      for (std::size_t i = 0; i < neighbours.size(); i++) {
        if (k + iatom * nlatticesites_ < neighbours[i])
          edge_vector.push_back({{k + iatom * nlatticesites_, neighbours[i]}});
      }
    }
  }
  return edge_vector;
};

std::vector<std::vector<int>> Lattice::AdjacencyList() const noexcept {
  std::vector<std::vector<int>> adlist(nsites_);
  for (int k = 0; k < nsites_; k++) {
    for (std::size_t i = 0; i < edges_.size(); i++) {
      if (edges_[i][0] == k)
        adlist[k].push_back(edges_[i][1]);
      else if (edges_[i][1] == k)
        adlist[k].push_back(edges_[i][0]);
    }
  }
  return adlist;
}

// Symmetries
std::vector<std::vector<int>> Lattice::BuildSymmetryTable() const {
  std::vector<std::vector<int>> integers_all;
  std::vector<std::vector<int>> result;
  for (int k = 0; k < nlatticesites_; k++)
    integers_all.push_back(Site2Vector(k));
  int i = 0;
  do {
    std::vector<int> row;
    int index;
    for (int iatom = 0; iatom < natoms_; iatom++) {
      for (int j = 0; j < nlatticesites_; j++) {
        std::vector<int> n(ndim_);
        index = 1;
        for (int l = 0; l < ndim_; l++) {
          if (pbc_[l]) {
            n[l] = (integers_all[i][l] + integers_all[j][l]) % (extent_[l]);
          } else {
            n[l] = integers_all[j][l];
            index *= extent_[l];
          }
        }
        row.push_back(Vector2Site(n) + iatom * nlatticesites_);  // save rows
      }
    }
    i += index;
    result.push_back(row);
  } while (i < nlatticesites_);

  return result;
};

bool Lattice::IsBipartite() const noexcept { return is_bipartite_; }

bool Lattice::IsConnected() const noexcept { return is_connected_; }

// Utilities
bool Lattice::RelativelyEqual(double a, double b,
                              double maxRelativeDiff) const {
  const double difference = std::abs(a - b);

  // Scale to the largest value.
  a = std::abs(a);
  b = std::abs(b);
  const double scaledEpsilon = maxRelativeDiff * std::max(a, b);
  return difference <= scaledEpsilon;
}
double Lattice::GetNorm(const std::vector<double> &coord) const {
  double distance = 0;
  for (std::size_t i = 0; i < coord.size(); i++)
    distance += coord[i] * coord[i];
  return distance;
};
double Lattice::GetSquaredDistance(const std::vector<double> &v1,
                                   const std::vector<double> &v2) const {
  double distance = 0;
  if (v1.size() != v2.size()) {
    ErrorMessage() << "Impossible to compute distance between two vectors of "
                      "different size. Exit.\n";
    exit(EXIT_FAILURE);
  }
  for (std::size_t i = 0; i < v1.size(); i++)
    distance += (v1[i] - v2[i]) * (v1[i] - v2[i]);
  return distance;
}

}  // namespace netket

#endif
