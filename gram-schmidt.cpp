#include "gram-schmidt.hpp"

// return the number of independent vectors in the basis
std::vector<Poly> Orthogonalize(const std::vector<Basis<Mono>>& inputBases, 
                std::ostream&) {
	Timer timer;
	Basis<Mono> unifiedBasis = CombineBases(inputBases);
	Normalize(unifiedBasis);
	// without the following stream, unifiedBasis segfaults
	// outStream << "Normalized initial basis: " << unifiedBasis << std::endl;

	DMatrix gram = GramFock(unifiedBasis);
	if(gram.rows() == 0) return {};
	
	std::cout << "Gram matrix constructed in " << timer.TimeElapsedInWords()
		<< "." << std::endl;
	if(TotalSize(inputBases) <= 7){
		std::cout << gram << std::endl;
	}

	// orthogonalize using custom gram-schmidt
	timer.Start();
	std::vector<Poly> orthogonalized = GramSchmidt_WithMatrix(unifiedBasis, gram);

	std::cout << "Gram-Schmidt performed in " << timer.TimeElapsedInWords()
		<< ", giving " << orthogonalized.size() << " vectors";
	if (orthogonalized.size() <= 20) {
		std::cout << ":" << std::endl;
		for(auto& p : orthogonalized) std::cout << p << std::endl;
	} else {
		std::cout << ", which will not be shown." << std::endl;
	}

    return orthogonalized;
}

std::vector<Poly> GramSchmidt_WithMatrix(const std::vector<Basis<Mono>> input,
		const DMatrix& gramMatrix) {
	return GramSchmidt_WithMatrix(CombineBases(input), gramMatrix);
}

/*std::vector<Poly> GramSchmidt_WithMatrix(const std::vector<Basis<Mono>> input,
		const DMatrixHP& gramMatrix) {
	return GramSchmidt_WithMatrix(CombineBases(input), gramMatrix);
}*/

// I've included a couple of implementations of the Gram-Schmidt algorithm; 
// these are equivalent for real numbers, but with floating point numbers 
// they'll produce different rounding errors. A seems to be better at the moment
// but neither one is entirely satisfactory.
std::vector<Poly> GramSchmidt_WithMatrix(const Basis<Mono> inputBasis, 
		const DMatrix& gramMatrix) {
	return GramSchmidt_WithMatrix_A(inputBasis, gramMatrix);
	// return GramSchmidt_WithMatrix_B (inputBasis, gramMatrix);
	// return GramSchmidt_WithMatrix_C (inputBasis, gramMatrix.cast<precise_class>());
}

/*std::vector<Poly> GramSchmidt_WithMatrix(const Basis<Mono> inputBasis, 
		const DMatrixHP& gramMatrix) {
	// return GramSchmidt_WithMatrix_A(inputBasis, gramMatrix.cast<coeff_class>());
	// return GramSchmidt_WithMatrix_B (inputBasis, gramMatrix.cast<coeff_class>());
	return GramSchmidt_WithMatrix_C (inputBasis, gramMatrix);
}*/

std::vector<Poly> GramSchmidt_WithMatrix_A(const Basis<Mono> inputBasis, 
		const DMatrix& gramMatrix) {
	std::vector<DVector> vectorForms;
	for (auto i = 0u; i < inputBasis.size(); ++i) {
		DVector nextVector = DVector::Unit(inputBasis.size(), i);
		for (auto j = 0u; j < vectorForms.size(); ++j) {
			nextVector -= GSProjection(nextVector, vectorForms[j], gramMatrix);
		}

		coeff_class norm = GSNorm(nextVector, gramMatrix);
		if (std::abs<builtin_class>(norm) < EPSILON) continue;
		if (norm < 0) {
			std::cerr << "Warning: negative norm " << norm << "." << std::endl;
			norm = -norm;
		}
		vectorForms.push_back(nextVector/std::abs(std::sqrt<builtin_class>(norm)));
	}

	std::vector<Poly> ret;
	for (auto& vec : vectorForms) ret.push_back(VectorToPoly(vec, inputBasis));
	return ret;
}

std::vector<Poly> GramSchmidt_WithMatrix_B(const Basis<Mono> inputBasis, 
		const DMatrix& gramMatrix) {
	std::vector<DVector> vectorForms;
	for (std::size_t i = 0; i < inputBasis.size(); ++i) {
		vectorForms.push_back(DVector::Unit(inputBasis.size(), i)
				/ std::abs(std::sqrt<builtin_class>(gramMatrix(i,i))) );
	}

	std::vector<coeff_class> norms(inputBasis.size(), 0);
	for (std::size_t i = 0; i < inputBasis.size(); ++i) {
		norms[i] = vectorForms[i].transpose() * gramMatrix * vectorForms[i];
		if (norms[i] < 0) norms[i] = -norms[i];
		if (norms[i] < EPSILON) continue;
		vectorForms[i] /= std::abs(std::sqrt<builtin_class>(norms[i]));
		for (std::size_t j = i+1; j < inputBasis.size(); ++j) {
			coeff_class projector =
				vectorForms[i].transpose() * gramMatrix * vectorForms[j];
			vectorForms[j] -= projector * vectorForms[i];
		}
	}

	std::vector<Poly> ret;
	for (std::size_t i = 0; i < inputBasis.size(); ++i) {
		if (norms[i] < EPSILON) continue;
		ret.push_back(VectorToPoly(vectorForms[i], inputBasis));
	}
	return ret;
}

// this is a clone of A but with higher precision for the coefficients
/*std::vector<Poly> GramSchmidt_WithMatrix_C(const Basis<Mono> inputBasis, 
		const DMatrixHP& gramMatrix) {
	std::vector<DVectorHP> vectorForms;
	for (auto i = 0u; i < inputBasis.size(); ++i) {
		DVectorHP nextVector = DVectorHP::Unit(inputBasis.size(), i);
		for (auto j = 0u; j < vectorForms.size(); ++j) {
			nextVector -= (nextVector.transpose() * gramMatrix * vectorForms[j])
				* vectorForms[j];
		}

		coeff_class norm = nextVector.transpose() * gramMatrix * nextVector;
		if (std::abs(norm) < EPSILON) continue;
		if (norm < 0) {
			std::cerr << "Warning: negative norm " << norm << "." << std::endl;
			norm = -norm;
		}
		vectorForms.push_back(nextVector/std::sqrt(norm));
	}

	std::vector<Poly> ret;
	for (auto& vec : vectorForms) ret.push_back(VectorToPoly(vec, inputBasis));
	return ret;
}*/

// projectOnto must be normalized already if you want the right answer from this
//
// the commented return line is what we're producing mathematically, but we
// can reduce rounding errors a bit by doing it by hand
DVector GSProjection(const DVector& toProject, const DVector& projectOnto,
		const DMatrix& gramMatrix) {
	coeff_class scale = 0;
	for (auto i = 0u; i < toProject.size(); ++i) {
		if (toProject(i) == 0) continue;
		for (auto j = 0u; j < toProject.size(); ++j) {
			scale += toProject(i)*projectOnto(j)*gramMatrix(i,j);
		}
	}
	return scale*projectOnto;
	// return (toProject.transpose() * gramMatrix * projectOnto) * projectOnto;
}

// the commented return line is what we're producing mathematically, but we
// can reduce rounding errors a bit by doing it by hand
coeff_class GSNorm(const DVector& vector, const DMatrix& gramMatrix) {
	coeff_class norm = 0;
	for(auto i = 0u; i < vector.size(); ++i) {
		norm += vector(i)*vector(i)*gramMatrix(i,i);
		for (auto j = i+1; j < vector.size(); ++j) {
			norm += 2*vector(i)*vector(j)*gramMatrix(i,j);
		}
	}
	return norm;
	// return vector.transpose() * gramMatrix * vector;
}

std::vector<Poly> GramSchmidt_MatrixOnly(const DMatrix& input, 
		const std::vector<Basis<Mono>>& inputBases) {
	std::vector<DVector> knownVectors;
	for (Eigen::Index col = 0; col < input.cols(); ++col) {
		DVector newVector = input.col(col);
		for (const DVector& knownVector : knownVectors) {
			newVector -= newVector.dot(knownVector)*knownVector;
		}
		if(newVector.dot(newVector) <= EPSILON) continue;
		newVector /= std::abs(std::sqrt<builtin_class>(newVector.dot(newVector)));
		knownVectors.push_back(newVector);
	}

	std::vector<Poly> output;
	for(const DVector& knownVector : knownVectors){
		Poly polyForm;
		//std::cout << "knownVector is an object with " << knownVector.rows()
			//<< " rows and " << knownVector.cols() << " columns." << std::endl;
		for(Eigen::Index i = 0; i < knownVector.rows(); ++i){
			polyForm += knownVector(i)*Get(inputBases, i);
		}
		output.push_back(polyForm);
	}
	return output;
}

// turns QMatrix into polynomials using using basis. gramMatrix is used to
// normalize the output, and rank is used to know how many to extract
std::vector<Poly> PolysFromQMatrix(const DMatrix& QMatrix, 
		const Basis<Mono>& basis, const DMatrix& gramMatrix,
		const Eigen::Index rank) {
	std::vector<Poly> output;
	for (Eigen::Index i = 0; i < rank; ++i) {
		Poly nextPoly;
		for (Eigen::Index j = 0; j < QMatrix.rows(); ++j) {
			nextPoly += QMatrix(j,i)*basis[j];
		}
		coeff_class norm = QMatrix.col(i).transpose()*gramMatrix*QMatrix.col(i);
		norm = std::abs(std::sqrt<builtin_class>(norm));
		output.push_back(nextPoly / norm);
	}
	// std::cout << QMatrix << "\nconverted to " << output << std::endl;
	return output;
}

DMatrix ExtractQMatrix(const Eigen::FullPivHouseholderQR<DMatrix>& solver, 
		               const int) {
	return solver.matrixQ();
}

DMatrix ExtractQMatrix(const Eigen::ColPivHouseholderQR<DMatrix>& solver, 
		               const int dimension) {
	return solver.householderQ()*DMatrix::Identity(dimension, dimension);
}
