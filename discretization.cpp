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

// take a matrix whose columns are the polynomials expressed in terms of a 
// minimal basis and discretize it
DMatrix DiscretizePolys(const DMatrix& polysOnMinBasis, 
        const std::size_t partitions) {
    DMatrix output = DMatrix::Zero(polysOnMinBasis.rows()*partitions,
                    polysOnMinBasis.cols()*partitions);
    for (Eigen::Index row = 0; row < polysOnMinBasis.rows(); ++row) {
        for (Eigen::Index col = 0; col < polysOnMinBasis.cols(); ++col) {
            for (std::size_t p = 0; p < partitions; ++p) {
                output(partitions*row + p, partitions*col + p) 
                    = polysOnMinBasis(row, col);
            }
        }
    }
    return output;
}
