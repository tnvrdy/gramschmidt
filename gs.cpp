/*
gs.cpp
------
we learned about the gram-schmidt process in linalg, and it's 
pretty algorithmic, so i decided to quickly implement it :)

input a set of nonzero n-vectors and get an orthogonal 
or orthonormal basis of the subspace they span!
*/

#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <numeric> // for inner_product (dot product)
#include <string>
#include <vector>
#include <cmath>

using namespace std;
using Vec = vector<double>;
using VecList = vector<vector<double>>;

bool haveSameDim(const VecList& vectors);
bool allNonzero(const VecList& vectors);
bool isNonzero(const Vec& vector);
void subtractEW(Vec& w, const Vec& projection);
Vec project(const Vec& g, const Vec& b);
void normalizeBasis(VecList& vectors);
void normalize(Vec& vector);
double magnitude(const Vec& vector);

// returns an orthogonal, or orthonormal, basis of the span of `vectors`. */
// `vectors` should contain 1 or more vectors, all nonzero, nonempty.
VecList gramSchmidt(const VecList& vectors, bool orthonormal) {
    if (!haveSameDim(vectors)) throw runtime_error("Input vectors must be of same dimension.");
    if (!allNonzero(vectors)) throw runtime_error("Input vectors must be nonzero.");

    VecList result = {vectors[0]};
    for (size_t i = 1; i < vectors.size(); i++) {
        const auto& v = vectors[i];
        auto w = v;
        for (size_t j = 0; j < i; j++) {
            subtractEW(w, project(v, result[j])); // classical GS! 
                                                  // modified GS would project wi, not vi, onto previous w's.
        }
        if (isNonzero(w)) result.push_back(w);
    }
    if (orthonormal) normalizeBasis(result);
    return result;
}

// returns whether or not all vectors in `vectors` have same dimension.
// assumes inputted vectors are nonempty.
bool haveSameDim(const VecList& vectors) {
    auto n = vectors[0].size();
    for (const auto& vector : vectors) {
        if (vector.size() != n) return false;
    }
    return true;
}

// returns whether or not all vectors in `vectors` are nonzero.
bool allNonzero(const VecList& vectors) {
    for (const auto& vector : vectors) {
        if (!isNonzero(vector)) return false;
    }
    return true;
}

// returns whether or not `vector` is nonzero.
bool isNonzero(const Vec& vector) {
    auto nonzero = false;
    for (double v : vector) {
        if (v != 0.0) {
            nonzero = true;
            break;
        }
    }
    return nonzero;
}

// element-wise subtraction of `p` from `w`.
void subtractEW(Vec& w, const Vec& p) {
    for (size_t i = 0; i < w.size(); i++) w[i] -= p[i];
}

// returns projection of `g` onto `b`.
Vec project(const Vec& g, const Vec& b) {
    Vec result(b.size());
    auto gDotB = inner_product(g.begin(), g.end(), b.begin(), 0.0);
    auto bDotB = inner_product(b.begin(), b.end(), b.begin(), 0.0);
    auto c = gDotB / bDotB;

    transform(b.begin(), b.end(), result.begin(),
              [c] (double v) {
                return v * c;
              });
    return result;
}

void normalizeBasis(VecList& vectors) {
    for (auto& vector : vectors) {
        normalize(vector);
    }
}

// transforms `vector` into its unit vector.
void normalize(Vec& vector) {
    auto m = magnitude(vector);
    transform(vector.begin(), vector.end(), vector.begin(),
              [m] (double v) {
                  return v / m;
              });
}

double magnitude(const Vec& vector) {
    auto sumSq = 0.0;
    for (double v : vector) sumSq += pow(v, 2);
    return sqrt(sumSq);
}

int main() { // leetol test
    VecList input = {{2.0, 1.0}, {1.0, 1.0}};
    auto result = gramSchmidt(input, true);

    for (size_t i = 0; i < result.size(); i++) {
        const auto& v = result[i];
        cout << "Vector " << i << ": { ";
        for (size_t j = 0; j < v.size(); j++) {
            cout << v[j] << " ";
        }
        cout << "}" << endl;
    }
}
