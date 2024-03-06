/*
 * mesh.h
 *   Finite element mesh implementation.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#ifndef MESH_H
#define MESH_H

#include "assembler.h"
#include "etype.h"

#include <vector>
using std::vector;

class Mesh {
public:
    Mesh(int ndm, int maxndf, int maxnen) :
	ndm_(ndm), maxndf_(maxndf), maxnen_(maxnen), numid_(0) {}
    ~Mesh();

    int ndm()    const { return ndm_;             }
    int maxndf() const { return maxndf_;          }
    int maxnen() const { return maxnen_;          }
    int numelt() const { return elements.size();  }
    int numnp()  const { return X.size() / ndm_;  }
    int numid()  const { return numid_;           }

    double* x () { return &(X[0]);  }
    int*    ix() { return &(IX[0]); }
    int*    id() { return &(ID[0]); }
    double* u () { return &(U[0]);  }
    double* f () { return &(R[0]);  }
    char*   bc() { return &(BC[0]); }
    double* bv() { return &(BV[0]); }

    const double* x (int j) const { return &(X [j*ndm_]);    }
    const int*    ix(int j) const { return &(IX[j*maxnen_]); }
    const int*    id(int j) const { return &(ID[j*maxndf_]); }
    const double* u (int j) const { return &(U [j*maxndf_]); }
    const double* f (int j) const { return &(R [j*maxndf_]); }
    double*       f (int j)       { return &(R [j*maxndf_]); }
    const char*   bc(int j) const { return &(BC[j*maxndf_]); }
    const double* bv(int j) const { return &(BV[j*maxndf_]); }

    double  x (int i, int j) const { return X [i+j*ndm_];    }
    int     ix(int i, int j) const { return IX[i+j*maxnen_]; }
    int     id(int i, int j) const { return ID[i+j*maxndf_]; }
    int&    id(int i, int j)       { return ID[i+j*maxndf_]; }
    double  u (int i, int j) const { return U [i+j*maxndf_]; }
    double  f (int i, int j) const { return R [i+j*maxndf_]; }
    double& f (int i, int j)       { return R [i+j*maxndf_]; }
    char    bc(int i, int j) const { return BC[i+j*maxndf_]; }
    char&   bc(int i, int j)       { return BC[i+j*maxndf_]; }
    double  bv(int i, int j) const { return BV[i+j*maxndf_]; }
    double& bv(int i, int j)       { return BV[i+j*maxndf_]; }

    void set_ur(const double* ur);
    void get_ur(double* ur);
    void get_fr(double* ur);

    int  initialize();
    int  assign_ids();

    void assemble_F();
    void assemble_K(MatrixAssembler* K_assembler);

    int  add_node(double* x);
    int  add_element(EType* etype, int* nodes, int numnp);
    void add_material(EType* material);

private:
    Mesh(const Mesh&);
    Mesh& operator=(const Mesh&);

    int ndm_;      // Number of spatial dimensions
    int maxndf_;   // Maximum degrees of freedom per node
    int maxnen_;   // Maximum number of nodes per element
    int numid_;    // Number of active degrees of freedom

    vector<double> X;   // Coordinate array      ( ndm    * numnp  )
    vector<int>    IX;  // Element connectivity  ( maxnen * numelt )
    vector<int>    ID;  // Identifier assignment ( maxndf * numnp  )
    vector<double> R;   // Full residual vector  ( maxndf * numnp  )
    vector<double> U;   // Full solution         ( maxndf * numnp  )
    vector<char>   BC;  // Flag Dirichlet BCs    ( maxndf * numnp  )
    vector<double> BV;  // Boundary values       ( maxndf * numnp  )

    vector<EType*> elements; // Material assignment (numelt)

    vector<EType*> materials_owned;
};

#endif /* MESH_H */
