#include "discretization.hpp"

// return exponent P of \mu^P to deal with when multiplying A and B; this is
// N-3 + total number of perp derivates + 2 because the integrals are originally
// over \mu^2 instead of \mu
char MuExponent(const Mono& A, const Mono& B) {
	return A.NParticles() - 3 + A.TotalPt() + B.TotalPt() + 2;
}

// this only returns a DVector instead of a DMatrix because the mass term has a
// \delta_{ij} so only the diagonal is nonzero for each pair of monomials
//
// the 4 is from converting two integrals from d\mu^2 to d\mu
DVector MuIntegral_Mass(const Mono& A, const Mono& B, 
		const std::size_t partitions, const coeff_class partitionWidth) {
	char muExp = MuExponent(A, B);
	DVector output(partitions);
	coeff_class left = 0;
	coeff_class right = 0;
	for (std::size_t i = 0; i < partitions; ++i) {
		left = right;
		right += partitionWidth;
		output[i] = 4 * (std::pow<double>(right, muExp+1) 
                        - std::pow<double>(left, muExp+1)) / (muExp + 1);
	}
	return output;
}

// this is a placeholder, obviously
DVector MuIntegral_InnerProduct(const Mono& A, const Mono& B, 
		const std::size_t partitions, const coeff_class partitionWidth) {
    return MuIntegral_Mass(A, B, partitions, partitionWidth);
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
			DVector integrated = MuIntegral_Mass(minimalBasis[row], 
					minimalBasis[col], partitions, partWidth);
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
    for (std::size_t i = 0; i < minBasis.size(); ++i) {
        DVector IP = MuIntegral_InnerProduct(
                    minBasis[i], minBasis[i], partitions, partWidth );
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

// take a matrix whose columns are the polynomials expressed in terms of a 
// minimal basis and discretize it
DMatrix DiscretizePolys(const DMatrix& polysOnMinBasis, 
        const Basis<Mono>& minBasis, const std::size_t partitions) {
    DMatrix output = DMatrix::Zero(polysOnMinBasis.rows()*partitions,
                    polysOnMinBasis.cols()*partitions);
    /*for (Eigen::Index row = 0; row < polysOnMinBasis.rows(); ++row) {
        for (Eigen::Index col = 0; col < polysOnMinBasis.cols(); ++col) {
            // this is minimalBasis[col] * minimalBasic[col] (no row) because
            // the row is indexing basis states, whereas the col is indexing
            // minBasis monomials
			DVector norms = MuIntegral_Mass(minimalBasis[col], 
					minimalBasis[col], partitions, partWidth);
            for (std::size_t p = 0; p < partitions; ++p) {
                output(partitions*row + p, partitions*col + p) 
                    = polysOnMinBasis(row, col) / std::sqrt(norms(p));
            }
        }
    }*/
    return output;
}
