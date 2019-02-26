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

#ifndef NETKET_PYLATTICE_HPP
#define NETKET_PYLATTICE_HPP

#include "lattice.hpp"

namespace py = pybind11;

namespace netket {

void AddLattice(py::module& subm) {
  py::class_<Lattice, AbstractGraph>(subm, "Lattice", R"EOF(
                             A generic lattice built translating a unit cell and adding edges between nearest neighbours sites.
                             The unit cell is defined by the ``basis_vectors`` and it can contain an arbitrary number of atoms.
                             Each atom is located at an arbitrary position and is labelled by an integer number,
                             meant to distinguish between the different atoms within the unit cell.
                             Periodic boundary conditions can also be imposed along the desired directions.
                             There are three different ways to refer to the lattice sites. A site can be labelled
                             by a simple integer number (the site index), by its coordinates (actual position in space),
                             or by a set of integers (the site vector), which indicates how many
                             translations of each basis vectors have been performed while building the
                             graph. The i-th component refers to translations along the i-th ``basis_vector`` direction.)EOF")
      .def(py::init<std::vector<std::vector<double>>, std::vector<int>,
                    std::vector<bool>, std::vector<std::vector<double>>>(),
           py::arg("basis_vectors"), py::arg("extent"),
           py::arg("pbc") = std::vector<bool>(0),
           py::arg("atoms_coord") = std::vector<std::vector<double>>(0),
           R"EOF(
                             Constructs a new ``Lattice`` given its side length and the features of the unit cell.

                             Args:
                                 basis_vectors: The basis vectors of the unit cell.
                                 extent: The number of copies of the unit cell.
                                 pbc: If ``True`` then the constructed lattice
                                     will have periodic boundary conditions, otherwise
                                     open boundary conditions are imposed (default=``True``).
                                 atoms_coord: The coordinates of different atoms in the unit cell (default=one atom at the origin).

                             Examples:
                                 Constructs a rectangular 3X4 lattice with periodic boundary conditions.

                                 ```python
                                 >>> from netket.graph import Lattice
                                 >>> g=Lattice(basis_vectors=[[1,0],[0,1]],extent=[3,4])
                                 >>> print(g.n_sites)
                                 12

                                 ```
                             )EOF")

      .def_property_readonly("coordinates", &Lattice::Coordinates,
                             R"EOF(
      list[list]: The coordinates of the atoms in the lattice.)EOF")
      .def_property_readonly("n_dim", &Lattice::Ndim, R"EOF(
      int: The dimension of the lattice.)EOF")
      .def_property_readonly("basis_vectors", &Lattice::BasisVectors, R"EOF(
      list[list]: The basis vectors of the lattice.)EOF")
      .def("atom_label", &Lattice::AtomLabel, py::arg("site"), R"EOF(
          Member function returning the atom label indicating which of the unit cell atoms is located at a given a site index.

          Args:
              site: The site index.

          )EOF")
      .def("site_to_vector", &Lattice::Site2Vector, py::arg("site"), R"EOF(
          Member function returning the site vector corresponding to a given site index.

          Args:
              site: The site index.

          Examples:
              Constructs a square 2X2 lattice without periodic boundary conditions and prints the site vectors corresponding to given site indices.

              ```python
               >>> from netket.graph import Lattice
               >>> g=Lattice(basis_vectors=[[1.,0.],[0.,1.]], extent=[2,2], pbc=[0,0])
               >>> print(list(map(int,g.site_to_vector(0))))
               [0, 0]
               >>> print(list(map(int,g.site_to_vector(1))))
               [0, 1]
               >>> print(list(map(int,g.site_to_vector(2))))
               [1, 0]
               >>> print(list(map(int,g.site_to_vector(3))))
               [1, 1]

                ```
            )EOF")
      .def("vector_to_coord", &Lattice::Vector2Coord, py::arg("site_vector"),
           py::arg("atom_label"), R"EOF(
        Member function returning the coordinates of a given atom characterized by a
        given site vector.

          Args:
              site_vector: The site vector.
              atom_label: Which of the atoms in the unit cell.

            )EOF")
      .def("site_to_coord", &Lattice::Site2Coord, py::arg("site"), R"EOF(
        Member function returning the coordinates of a given site index.

          Args:
              site: The site index.
            )EOF")
      .def("vector_to_site", &Lattice::Vector2Site, py::arg("site_vector"),
           R"EOF(
        Member function returning the site index corresponding to a given site vector.

            Args:
                site_vector: The site vector.
            )EOF")
      // Built-in lattices
      .def_static("hypercube",
                  [](int dim, std::vector<int> extent,
                     std::vector<bool> pbc) -> Lattice {
                    std::vector<std::vector<double>> basis_vectors(
                        dim, std::vector<double>(dim));
                    std::vector<std::vector<double>> atoms_coord(
                        1, std::vector<double>(dim));
                    for (int i = 0; i < dim; i++) {
                      atoms_coord[0][i] = 0;
                      for (int j = 0; j < dim; j++) {
                        if (i == j) {
                          basis_vectors[i][j] = 1;
                        } else {
                          basis_vectors[i][j] = 0;
                        }
                      }
                    }
                    return Lattice(basis_vectors, extent, pbc, atoms_coord);
                  },
                  py::arg("n_dim"), py::arg("extent"),
                  py::arg("pbc") = std::vector<bool>(0), R"EOF(
                Member function constructing a hypercubic lattice of arbitrary
                dimension.

                  Args:
                      n_dim: The dimension of the lattice.
                      extent: The number of copies of the unit cell.
                      pbc: If ``True`` then the constructed lattice will have periodic boundary conditions, otherwise open boundary conditions are imposed (default=``True``).
                  )EOF")

      .def_static(
          "monoclinic",
          [](std::vector<int> extent, std::vector<bool> pbc) -> Lattice {
            std::vector<std::vector<double>> basis_vectors(
                2, std::vector<double>(2));
            std::vector<std::vector<double>> atoms_coord(
                1, std::vector<double>(2));
            basis_vectors[0][0] = 1.;
            basis_vectors[0][1] = 0.;
            basis_vectors[1][0] = 1. / 3.;
            basis_vectors[1][1] = 2. / 3.;
            atoms_coord[0][0] = 0;
            atoms_coord[0][1] = 0;
            return Lattice(basis_vectors, extent, pbc, atoms_coord);
          },
          py::arg("extent"), py::arg("pbc") = std::vector<bool>(0), R"EOF(
      Member function constructing a monoclinic lattice in 2 dimensions.

        Args:
            extent: The number of copies of the unit cell.
            pbc: If ``True`` then the constructed lattice will have periodic boundary conditions, otherwise open boundary conditions are imposed (default=``True``).


        )EOF")
      .def_static(
          "rectangular",
          [](std::vector<int> extent, std::vector<bool> pbc) -> Lattice {
            std::vector<std::vector<double>> basis_vectors(
                2, std::vector<double>(2));
            std::vector<std::vector<double>> atoms_coord(
                1, std::vector<double>(2));
            basis_vectors[0][0] = 1.;
            basis_vectors[0][1] = 0.;
            basis_vectors[1][0] = 0.;
            basis_vectors[1][1] = 2.;
            atoms_coord[0][0] = 0;
            atoms_coord[0][1] = 0;

            return Lattice(basis_vectors, extent, pbc, atoms_coord);
          },
          py::arg("extent"), py::arg("pbc") = std::vector<bool>(0), R"EOF(
      Member function constructing a rectangular lattice in 2 dimensions.

        Args:
            extent: The number of copies of the unit cell.
            pbc: If ``True`` then the constructed lattice will have periodic boundary conditions, otherwise open boundary conditions are imposed (default=``True``).


        )EOF")

      .def_static(
          "centered_rectangular",
          [](std::vector<int> extent, std::vector<bool> pbc) -> Lattice {
            std::vector<std::vector<double>> basis_vectors(
                2, std::vector<double>(2));
            std::vector<std::vector<double>> atoms_coord(
                2, std::vector<double>(2));
            basis_vectors[0][0] = 1.;
            basis_vectors[0][1] = 0.;
            basis_vectors[1][0] = 0.;
            basis_vectors[1][1] = 2.;
            atoms_coord[0][0] = 0.;
            atoms_coord[0][1] = 0.;
            atoms_coord[1][0] = 0.5;
            atoms_coord[1][1] = 1.;

            return Lattice(basis_vectors, extent, pbc, atoms_coord);
          },
          py::arg("extent"), py::arg("pbc") = std::vector<bool>(0), R"EOF(
      Member function constructing a centered rectangular lattice in 2 dimensions.

        Args:
            extent: The number of copies of the unit cell.
            pbc: If ``True`` then the constructed lattice will have periodic boundary conditions, otherwise open boundary conditions are imposed (default=``True``).


        )EOF")

      .def_static("hexagonal",
                  [](int dim, std::vector<int> extent,
                     std::vector<bool> pbc) -> Lattice {
                    std::vector<std::vector<double>> basis_vectors(
                        dim, std::vector<double>(dim));
                    std::vector<std::vector<double>> atoms_coord(
                        1, std::vector<double>(dim));
                    if (dim == 2) {
                      basis_vectors[0][0] = 1.;
                      basis_vectors[0][1] = 0.;
                      basis_vectors[1][0] = 0.5;
                      basis_vectors[1][1] = 0.5 * std::sqrt(3);
                      atoms_coord[0][0] = 0;
                      atoms_coord[0][1] = 0;
                    } else if (dim == 3) {
                      basis_vectors[0][0] = 1.;
                      basis_vectors[0][1] = 0.;
                      basis_vectors[0][2] = 0.;
                      basis_vectors[1][0] = 0.5;
                      basis_vectors[1][1] = 0.5 * std::sqrt(3);
                      basis_vectors[1][2] = 0.;
                      basis_vectors[2][0] = 0.;
                      basis_vectors[2][1] = 0.;
                      basis_vectors[2][2] = 1.;
                      atoms_coord[0][0] = 0;
                      atoms_coord[0][1] = 0;
                      atoms_coord[0][2] = 0;
                    } else {
                      throw InvalidInputError{
                          "Lattice.hexagonal initializer works only in 2 and 3 "
                          "dimensions.\n"};
                    }
                    return Lattice(basis_vectors, extent, pbc, atoms_coord);
                  },
                  py::arg("n_dim"), py::arg("extent"),
                  py::arg("pbc") = std::vector<bool>(0), R"EOF(
    Member function constructing a hexagonal lattice in 2 or 3 dimensions.

      Args:
          n_dim: The dimension of the lattice. It can be 2 or 3.
          extent: The number of copies of the unit cell.
          pbc: If ``True`` then the constructed lattice will have periodic boundary conditions, otherwise open boundary conditions are imposed (default=``True``).


      )EOF")

      .def_static(
          "honeycomb",
          [](std::vector<int> extent, std::vector<bool> pbc) -> Lattice {
            std::vector<std::vector<double>> basis_vectors(
                2, std::vector<double>(2));
            std::vector<std::vector<double>> atoms_coord(
                2, std::vector<double>(2));
            basis_vectors[0][0] = 1.5;
            basis_vectors[0][1] = std::sqrt(3) * 0.5;
            basis_vectors[1][0] = 0.;
            basis_vectors[1][1] = std::sqrt(3);
            atoms_coord[0][0] = 0.;
            atoms_coord[0][1] = 0.;
            atoms_coord[1][0] = 1.;
            atoms_coord[1][1] = 0.;

            return Lattice(basis_vectors, extent, pbc, atoms_coord);
          },
          py::arg("extent"), py::arg("pbc") = std::vector<bool>(0), R"EOF(
    Member function constructing a honeycomb lattice in 2 dimensions.

      Args:
          extent: The number of copies of the unit cell.
          pbc: If ``True`` then the constructed lattice will have periodic boundary conditions, otherwise open boundary conditions are imposed (default=``True``).


      )EOF")

      .def_static(
          "kagome",
          [](std::vector<int> extent, std::vector<bool> pbc) -> Lattice {
            std::vector<std::vector<double>> basis_vectors(
                2, std::vector<double>(2));
            std::vector<std::vector<double>> atoms_coord(
                3, std::vector<double>(2));

            basis_vectors[0][0] = 2.;
            basis_vectors[0][1] = 0.;
            basis_vectors[1][0] = 1.;
            basis_vectors[1][1] = std::sqrt(3);
            atoms_coord[0][0] = 0.;
            atoms_coord[0][1] = 0.;
            atoms_coord[1][0] = 0.5;
            atoms_coord[1][1] = std::sqrt(3) * 0.5;
            atoms_coord[2][0] = 1.;
            atoms_coord[2][1] = 0.;

            return Lattice(basis_vectors, extent, pbc, atoms_coord);
          },
          py::arg("extent"), py::arg("pbc") = std::vector<bool>(0), R"EOF(
  Member function constructing a kagome lattice in 2 dimensions.

    Args:
        extent: The number of copies of the unit cell.
        pbc: If ``True`` then the constructed lattice will have periodic boundary conditions, otherwise open boundary conditions are imposed (default=``True``).


    )EOF")
      .def_static(
          "triclinic",
          [](std::vector<int> extent, std::vector<bool> pbc) -> Lattice {
            std::vector<std::vector<double>> basis_vectors(
                3, std::vector<double>(3));
            std::vector<std::vector<double>> atoms_coord(
                1, std::vector<double>(3));

            basis_vectors[0][0] = 1.;
            basis_vectors[0][1] = 0.;
            basis_vectors[0][2] = 0.;
            basis_vectors[1][0] = -0.5;
            basis_vectors[1][1] = 0.;
            basis_vectors[1][2] = 0.5;
            basis_vectors[2][0] = 1.;
            basis_vectors[2][1] = 1.;
            basis_vectors[2][2] = 0.;
            atoms_coord[0][0] = 0.;
            atoms_coord[0][1] = 0.;
            atoms_coord[0][2] = 0.;

            return Lattice(basis_vectors, extent, pbc, atoms_coord);
          },
          py::arg("extent"), py::arg("pbc") = std::vector<bool>(0), R"EOF(
      Member function constructing a triclinic lattice in 3 dimensions.

        Args:
            extent: The number of copies of the unit cell.
            pbc: If ``True`` then the constructed lattice will have periodic boundary conditions, otherwise open boundary conditions are imposed (default=``True``).


        )EOF")

      .def_static("orthorhombic",
                  [](std::vector<int> extent, std::vector<bool> pbc,
                     std::string type) -> Lattice {
                    std::vector<std::vector<double>> basis_vectors(
                        3, std::vector<double>(3));
                    int natoms;
                    if (type.empty()) type == "primitive";
                    if (type == "primitive") {
                      natoms = 1;
                    } else if (type == "base") {
                      natoms = 2;
                    } else if (type == "bc") {
                      natoms = 2;
                    } else if (type == "fc") {
                      natoms = 4;
                    }
                    std::vector<std::vector<double>> atoms_coord(
                        natoms, std::vector<double>(3));

                    basis_vectors[0][0] = 1.;
                    basis_vectors[0][1] = 0.;
                    basis_vectors[0][2] = 0.;
                    basis_vectors[1][0] = 0.;
                    basis_vectors[1][1] = 2.;
                    basis_vectors[1][2] = 0.;
                    basis_vectors[2][0] = 0.;
                    basis_vectors[2][1] = 0.;
                    basis_vectors[2][2] = 0.5;
                    if (type == "primitive") {
                      atoms_coord[0][0] = 0;
                      atoms_coord[0][1] = 0;
                      atoms_coord[0][2] = 0;
                    } else if (type == "base") {
                      atoms_coord[0][0] = 0;
                      atoms_coord[0][1] = 0;
                      atoms_coord[0][2] = 0;
                      atoms_coord[1][0] = 0.5;
                      atoms_coord[1][1] = 1.;
                      atoms_coord[1][2] = 0.;
                    } else if (type == "bc") {
                      atoms_coord[0][0] = 0;
                      atoms_coord[0][1] = 0;
                      atoms_coord[0][2] = 0;
                      atoms_coord[1][0] = 0.5;
                      atoms_coord[1][1] = 1.;
                      atoms_coord[1][2] = 0.25;
                    } else if (type == "fc") {
                      atoms_coord[0][0] = 0;
                      atoms_coord[0][1] = 0;
                      atoms_coord[0][2] = 0;
                      atoms_coord[1][0] = 0.5;
                      atoms_coord[1][1] = 1.;
                      atoms_coord[1][2] = 0.;
                      atoms_coord[2][0] = 0.5;
                      atoms_coord[2][1] = 0.;
                      atoms_coord[2][2] = 0.25;
                      atoms_coord[3][0] = 0.;
                      atoms_coord[3][1] = 1.;
                      atoms_coord[3][2] = 0.25;
                    }

                    return Lattice(basis_vectors, extent, pbc, atoms_coord);
                  },
                  py::arg("extent"), py::arg("pbc") = std::vector<bool>(0),
                  py::arg("type") = "primitive", R"EOF(
    Member function constructing an orthorhombic lattice in 3 dimensions.

      Args:
          extent: The number of copies of the unit cell.
          pbc: If ``True`` then the constructed lattice will have periodic boundary conditions, otherwise open boundary conditions are imposed (default=``True``).
          type: Wheter the lattice is primitive (default), base-centered (``base``), body-centered (``bc``) or face-centered (``fc``).

      )EOF");
}
}  // namespace netket
#endif
