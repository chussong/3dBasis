#include "discretization.hpp"

// return exponent P of \mu^P to deal with when multiplying A and B; this is
// N-3 + total number of perp derivates + 1 because the integrals are originally
// over \mu^2 instead of \mu (and I think only one of them survives so not +2??)
//
// I believe this corresponds to the exponent of \mu found in the inner product
char MuExponent(const Mono& A, const Mono& B) {
	return A.NParticles() - 3 + A.TotalPt() + B.TotalPt() + 1;
}

DVector MuIntegral(const Mono& A, const Mono& B, const std::size_t partitions,
        const coeff_class partWidth, const MATRIX_TYPE calculationType) {
    if (calculationType == MAT_INNER || calculationType == MAT_MASS) {
        DVector output = MuIntegral_Body(MuExponent(A, B), partitions, partWidth); 
        for (Eigen::Index k = 0; k < output.size(); ++k) {
            output(k) *= MuNorm(A, k, partWidth)*MuNorm(B, k, partWidth)/M_PI;
        }
        return output;
    } else {
        return DVector::Zero(partitions);
    }
}

// this only returns a DVector instead of a DMatrix because the mass term has a
// \delta_{ij} so only the diagonal is nonzero for each pair of monomials
//
// the 4 is from converting two integrals from d\mu^2 to d\mu
DVector MuIntegral_Body(const char muExp, const std::size_t partitions,
                const coeff_class partitionWidth) {
    DVector output(partitions);
    coeff_class left = 0;
    coeff_class right = 0;
    for (std::size_t i = 0; i < partitions; ++i) {
        left = right;
        right += partitionWidth;
        output[i] = (std::pow<double>(right, muExp+1) 
                - std::pow<double>(left, muExp+1)) / (muExp + 1);
        // std::cout << "Integral(" << (int)muExp << ", " << i << ", " 
            // << partitionWidth << ") = " << output[i] << std::endl;
    }
    return output;
}

// normalization for mu integrals; this is the script N_{n, \lambda, k} in
// Zuhair's notes
coeff_class MuNorm(const Mono& A, const std::size_t k, 
        const coeff_class partWidth) {
    // auto exponent = A.NParticles() - 1 + 2*A.TotalPt();
    auto exponent = MuExponent(A, A) + 1;
    coeff_class left = k*partWidth;
    coeff_class right = (k+1)*partWidth;

    // coeff_class denom = std::pow<builtin_class>(right, exponent);
    // denom -= std::pow<builtin_class>(left, exponent);
    // denom /= M_PI*exponent;

    coeff_class square = M_PI*exponent/(std::pow<builtin_class>(right, exponent)
            - std::pow<builtin_class>(left, exponent));

    auto output = std::sqrt<builtin_class>(square);
    if (output.imag() != 0) {
        std::cerr << "Warning: MuNorm(" << A << ", " << k << ", " << partWidth
            << ") is complex!?" << std::endl;
    }
    // std::cout << "MuNorm(" << A << ", " << k << ", " << partWidth << ") = "
        // << output.real() << std::endl;
    return output.real()/2;
}

// expands a mass matrix with the mu factors elided into one with them partially
// integrated over, by taking each entry of the input matrix and turning it into
// a block. Each block is diagonal because this is a mass matrix; it's not quite
// so simple for the interaction term
DMatrix PartitionMu_Mass(const Basis<Mono>& minimalBasis, const DMatrix& mass,
		const std::size_t partitions, const coeff_class partWidth) {
    DMatrix output = DMatrix::Zero(mass.rows()*partitions, mass.cols()*partitions);
    for (Eigen::Index row = 0; row < mass.rows(); ++row) {
        for (Eigen::Index col = 0; col < mass.cols(); ++col) {
            DVector integrated = MuIntegral(minimalBasis[row], 
                            minimalBasis[col], partitions, partWidth,
                            MAT_MASS);
            for (std::size_t p = 0; p < partitions; ++p) {
                    output(row*partitions + p, col*partitions + p) 
                            = mass(row, col) * integrated(p);
            }
        }
    }
    return output;
}

// for each monomial in the basis, create a vector with 1/sqrt(IPintegral) over
// each mu window. Return these as a matrix; the full basis state will then be
// [\sum_b X_{ib} Y_{bp} M_b], where X is the polysOnMinBasis matrix, Y is the
// output of this function, and M is the minBasis
DMatrix DiscretizeMonos(const Basis<Mono>& minBasis, 
        const std::size_t partitions, const coeff_class partWidth) {
    DMatrix output(minBasis.size(), partitions);
    /**************************************************************************/
    /* DANGER! PRETTY SURE THIS SHOULD BE DELETED AND IS NOT USEFUL AT ALL!!! */
    /**************************************************************************/
    for (std::size_t i = 0; i < minBasis.size(); ++i) {
        DVector IP = MuIntegral(
                    minBasis[i], minBasis[i], partitions, partWidth, MAT_INNER);
        for (std::size_t p = 0; p < partitions; ++p) {
            auto root = std::sqrt<builtin_class>(IP(p));
            if (root.imag() != 0) {
                throw std::logic_error("DiscretizeMonos: negative IP given.");
            }
            output(i, p) = 1/root.real();
        }
    }
    return output;
}

// same as above, except instead of a rank 2 tensor giving "integral over
// partition p" it's a rank 3 tensor giving "integral with monomial b over
// partition p"
DMatrix MuPart(const Basis<Mono>& minBasis, const std::size_t partitions, 
        const coeff_class partWidth, MATRIX_TYPE calculationType) {
    // expand discretized mono matrix into one that transforms polysOnMinBasis
    // (rank 2) into polysOnDiscretizedMinBasis, the rank 3 tensor above. This
    // means translating "partition p" into "partition p of monomial b".
    DMatrix output(minBasis.size(), minBasis.size()*partitions);
    for (std::size_t monoA = 0; monoA < minBasis.size(); ++monoA) {
        for (std::size_t monoB = 0; monoB < minBasis.size(); ++monoB) {
            DVector integral = MuIntegral(minBasis[monoA], minBasis[monoB], 
                    partitions, partWidth, calculationType);
            for (std::size_t part = 0; part < partitions; ++part) {
                output(monoA, monoB*partitions + part) = integral(part);
            }
        }
    }
    return output;
}

// return the mu integral, summed over all partitions, for each pair of
// monomials in the given basis
DMatrix MuTotal(const Basis<Mono>& minBasis, const std::size_t partitions,
        const coeff_class partWidth, const MATRIX_TYPE calculationType) {
    DMatrix output(minBasis.size(), minBasis.size());
    for (Eigen::Index row = 0; row < output.rows(); ++row) {
        DVector windows = MuIntegral(minBasis[row], minBasis[row], partitions,
                partWidth, calculationType);
        output(row, row) = windows.sum();
        for (Eigen::Index col = row+1; col < output.cols(); ++col) {
            DVector windows = MuIntegral(minBasis[row], minBasis[col], 
                    partitions, partWidth, calculationType);
            output(row, col) = windows.sum();
            output(col, row) = output(row, col);
        }
    }
    return output;
}

// take a matrix whose columns are the polynomials expressed in terms of a 
// minimal basis and discretize it according to the above procedure.
//
// each basis state is S_{i,p} = [\sum_b X_{ib} Y_{bp} M_{bp}], where X is 
// polysOnMinBasis, Y is the matrix from DiscretizeMonos (discMinBasis), 
// and M is minBasis split into partitions
//
// the matrix returned from this is effectively the rank 3 tensor S_{ibp},
// because the column index corresponds to j = P*b + p, where P is the total
// number of partitions. This means we should contract it with a rank 2 tensor
// M_{bp} in order to get a vector of S_i
DMatrix DiscretizePolys(const DMatrix& polysOnMinBasis, 
        const Basis<Mono>& minBasis, const std::size_t partitions,
        const coeff_class partWidth) {
    // expand discretized mono matrix into one that transforms polysOnMinBasis
    // (rank 2) into polysOnDiscretizedMinBasis, the rank 3 tensor above. This
    // means translating "partition p" into "partition p of monomial b".
    //
    // in other words: if this matrix is Y, Y_{bb'p} = MuIntegral(b, b', p)

    DMatrix matrixY = DMatrix::Zero(minBasis.size(), minBasis.size()*partitions);
    for (std::size_t i = 0; i < minBasis.size(); ++i) {
        // for (std::size_t j = 0; j < minBasis.size(); ++j) {
            // DVector integrals = MuIntegral(minBasis[i], minBasis[j], 
                    // partitions, partWidth, MAT_INNER);
            // for (std::size_t part = 0; part < partitions; ++part) {
                // matrixY(i, j*partitions + part) = integrals(part);
            // }
        // }
        for (std::size_t part = 0; part < partitions; ++part) {
            matrixY(i, i*partitions + part) = 1;
        }
    }

    return polysOnMinBasis * matrixY;
}
