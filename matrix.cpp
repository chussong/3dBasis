#include "matrix.hpp"

// Fock space inner product between two monomials
coeff_class InnerFock(const Mono& A, const Mono& B) {
	//B.ChangePm(0, B.Pm(0)+1); // this will be reverted in MatrixTerm below
	return MatrixInternal::MatrixTerm(A, B, MatrixInternal::MAT_INNER);
	//return MatrixInternal::MatrixTerm(A, B, MatrixInter::MAT_INNER)/A.NParticles();
}

// creates a gram matrix for the given basis using the Fock space inner product
DMatrix GramFock(const Basis<Mono>& basis) {
	DMatrix gramMatrix(basis.size(), basis.size());
	for (std::size_t i = 0; i < basis.size(); ++i) {
		gramMatrix(i, i) = InnerFock(basis[i], basis[i]);
		for (std::size_t j = i+1; j < basis.size(); ++j) {
			gramMatrix(i, j) = InnerFock(basis[i], basis[j]);
			gramMatrix(j, i) = gramMatrix(i, j);
		}
	}
	return gramMatrix;
}

// creates a mass matrix M for the given monomials. To get the mass matrix of a 
// basis of primary operators, one must express the primaries as a matrix of 
// vectors, A, and multiply A^T M A.
DMatrix MassMatrix(const Basis<Mono>& basis) {
	DMatrix massMatrix(basis.size(), basis.size());
	for (std::size_t i = 0; i < basis.size(); ++i) {
		massMatrix(i, i) = MatrixInternal::MatrixTerm(basis[i], basis[i], 
													MatrixInternal::MAT_MASS);
		for (std::size_t j = i+1; j < basis.size(); ++j) {
			massMatrix(i, j) = MatrixInternal::MatrixTerm(basis[i], basis[j], 
													MatrixInternal::MAT_MASS);
			massMatrix(j, i) = massMatrix(i, j);
		}
	}
	return massMatrix;
}

namespace MatrixInternal {

YTerm::YTerm(const coeff_class coeff, const std::string& y, 
		const std::string& nAndm): coeff(coeff), y(y.begin(), y.end()-1) {
	for (std::size_t i = 0; i < size(); ++i) {
		this->y[i] += nAndm[i+1];
	}
}

std::ostream& operator<<(std::ostream& os, const YTerm& out) {
	return os << out.coeff << " * " << out.y;
}

MatrixTerm_Intermediate::MatrixTerm_Intermediate(const size_t n): coefficient(1),
	uPlus(n), uMinus(n), yTilde(n) {
}

void MatrixTerm_Intermediate::Resize(const size_t n) {
	uPlus.resize(n);
	uMinus.resize(n);
	yTilde.resize(n);
}

MatrixTerm_Intermediate operator*(
		const MatrixTerm_Intermediate& A, MatrixTerm_Intermediate B) {
	B.coefficient *= A.coefficient;
	B.uPlus  = AddVectors(A.uPlus , B.uPlus );
	B.uMinus = AddVectors(A.uMinus, B.uMinus);
	B.yTilde = AddVectors(A.yTilde, B.yTilde);
	return B;
}

std::ostream& operator<<(std::ostream& os, const MatrixTerm_Intermediate& out) {
	return os << out.coefficient << " * {" << out.uPlus << ", "
		<< out.uMinus << ", " << out.yTilde << "}";
}

// this n would be called (n-1) in Zuhair's equations
MatrixTerm_Final::MatrixTerm_Final(const size_t n): coefficient(1), uPlus(n), 
	uMinus(n), sinTheta(n-1), cosTheta(n-1) {
}

// maybe should be rvalue references instead? hopefully it's the same
MatrixTerm_Final::MatrixTerm_Final(const coeff_class coefficient, 
		const std::vector<char>& uPlus, const std::vector<char>& uMinus, 
		const std::vector<char>& sinTheta, const std::vector<char>& cosTheta): 
	coefficient(coefficient), uPlus(uPlus), uMinus(uMinus), sinTheta(sinTheta), 
	cosTheta(cosTheta) {
}

// notice that there are fewer Thetas than Us, so the vectors aren't all the
// same size; also, the n here would be called (n-1) in Zuhair's equations
void MatrixTerm_Final::Resize(const size_t n) {
	uPlus.resize(n);
	uMinus.resize(n);
	sinTheta.resize(n-1);
	cosTheta.resize(n-1);
}

coeff_class MatrixTerm(const Mono& A, const Mono& B, 
		const MATRIX_TYPE type) {
	//std::cout << "TERM: " << A.HumanReadable() << " x " << B.HumanReadable() 
		//<< std::endl;
	//B.ChangePm(0, B.Pm(0)-1);

	// degeneracy factors result from turning the ordered monomials into 
	// symmetric polynomials
	coeff_class degeneracy = 1;
	// I think the A step here is wrong and I have to instead be weighting each
	// total configuration by its A degeneracy?
	// for (auto& count : A.CountIdentical()) degeneracy *= Factorial(count);
	degeneracy *= Factorial(A.NParticles());
	for (auto& count : B.CountIdentical()) degeneracy *= Factorial(count);

	std::vector<std::size_t> permA(A.PermutationVector());
	std::vector<std::size_t> permB(B.PermutationVector());
	std::array<std::string,2> xAndy_A, xAndy_B;
	std::vector<MatrixTerm_Final> fFromA, fFromB, combinedFs;
	coeff_class total = 0;
	// there's no reason to be using these permutation vectors instead of 
	// permuting xAndy directly
	//
	// one of these permutations is likely unnecesarry, which should be a
	// significant time savings
	// do {
		fFromA = MatrixTermsFromMono_Permuted(A, permA);
		// xAndy_A = ExtractXY(A);
		// xAndy_A = PermuteXandY(xAndy_A, permA);
		do {
			fFromB = MatrixTermsFromMono_Permuted(B, permB);
			combinedFs = CombineTwoFs(fFromA, fFromB);
			// std::array<std::string,2> xAndy_B = ExtractXY(B);
			// xAndy_B = PermuteXandY(xAndy_B, permB);
			// xAndy_B = CombineXandY(xAndy_A, xAndy_B);
			// combinedFs = MatrixTermsFromXandY(xAndy_B, A.NParticles());
			total += FinalResult(combinedFs, type);
		} while (std::next_permutation(permB.begin(), permB.end()));
	// } while (std::next_permutation(permA.begin(), permA.end()));
	
	if (type == MAT_INNER) {
		total *= InnerProductPrefactor(A.NParticles());
	} else if (type == MAT_MASS) {
		total *= MassMatrixPrefactor(A.NParticles());
	} else if (type == MAT_PHI4) {
	}

	// this somewhat dubious adjustment is presumably due to an error in
	// Zuhair's formula (the adjustment appears in his code as well)
	if (A.NParticles() >= 3) total /= 2;

	return degeneracy*A.Coeff()*B.Coeff()*total;
}

// takes a Mono and gets the exponents of the u and theta variables times some
// coefficient; the vector returned here contains the components of the sum
// constituting the answer, each with its own exponents and coefficient
std::vector<MatrixTerm_Final> MatrixTermsFromMono(const Mono& input) {
	std::array<std::string,2> xAndy(ExtractXY(input));
	return MatrixTermsFromXandY(xAndy, input.NParticles());
}

// this is a silly function and we should just be generating xAndy outside
// and permuting it directly instead of calling this over and over
std::vector<MatrixTerm_Final> MatrixTermsFromMono_Permuted(const Mono& input,
		const std::vector<std::size_t>& permutationVector) {
	std::array<std::string,2> xAndy(ExtractXY(input));
	//std::cout << "Extracted " << xAndy[0] << " from " << input;
	xAndy = PermuteXandY(xAndy, permutationVector);
	//std::cout << " and permuted it to " << xAndy[0] << std::endl;
	return MatrixTermsFromXandY(xAndy, input.NParticles());
}

// !!! NOTE: THE RESULTS OF THIS ENTIRE FUNCTION SHOULD BE HASHED !!!
std::vector<MatrixTerm_Final> MatrixTermsFromXandY(
		const std::array<std::string,2>& xAndy, const int nParticles) {
	std::vector<char> uFromX(UFromX(xAndy[0]));
	//std::cout << "uFromX: " << uFromX << std::endl;
	// !!! NOTE: SHOULD HASH THE RESULTS OF BELOW FUNCTIONS !!!
	// std::vector<MatrixTerm_Intermediate> inter(YTildeFromY(xAndy[1]));
	// std::vector<MatrixTerm_Final> ret(ThetaFromYTilde(inter));
	std::vector<MatrixTerm_Final> ret(ThetaFromY(xAndy[1]));
	// u has contributions from both x and y, so we have to combine them
	if (ret.size() == 0) ret.emplace_back(nParticles - 1);
	for (auto& term : ret) {
		if (term.uPlus.size() < uFromX.size()/2) {
			term.uPlus.resize(uFromX.size()/2, 0);
			term.uMinus.resize(uFromX.size()/2, 0);
			term.sinTheta.resize(uFromX.size()/2 - 1, 0);
			term.cosTheta.resize(uFromX.size()/2 - 1, 0);
		}
		for (auto i = 0u; i < term.uPlus.size(); ++i) {
			term.uPlus[i] += uFromX[i];
			term.uMinus[i] += uFromX[term.uPlus.size() + i];
		}
	}
	if (ret.size() == 0) {
		std::cerr << "Warning: returning an empty set of MatrixTermsFromMono. "
			<< "an empty set should contain one with coefficient 1 and all "
			<< "exponents set to 0." << std::endl;
	}
	return ret;
}

// exponent transformations ---------------------------------------------------
//
// rather than transforming the variables themselves, these take in a list of 
// exponents in one coordinate system and give out lists of exponents in the new 
// coordinate system.

// just gets the exponents of each x and y, so we don't have to worry about
// things like dividing by P or \mu
//
// note that I'm storing the Dirichlet-mandated P_- on each particle but Zuhair
// isn't, so we have to subtract that off
std::array<std::string,2> ExtractXY(const Mono& extractFromThis) {
	std::string x, y;
	for (auto i = 0u; i < extractFromThis.NParticles(); ++i) {
		x.push_back(extractFromThis.Pm(i) - 1);
		y.push_back(extractFromThis.Pt(i));
	}
	return {{x, y}};
}

// this function is pointless; xAndy should be generated once and permuted
// directly instead of remaking it over and over again
std::array<std::string,2> PermuteXandY(
		const std::array<std::string,2>& xAndy,
		const std::vector<std::size_t>& permutationVector) {
	std::array<std::string,2> output;
	output[0].resize(xAndy[0].size());
	output[1].resize(xAndy[1].size());
	for (std::size_t i = 0; i < permutationVector.size(); ++i) {
		output[0][i] = xAndy[0][permutationVector[i]];
		output[1][i] = xAndy[1][permutationVector[i]];
	}
	//std::cout << "Permuted " << xAndy[0] << " to " << output[0] << " using "
		//<< permutationVector << std::endl;
	return output;
}

std::array<std::string,2> CombineXandY(const std::array<std::string,2>& xAndy_A,
		std::array<std::string,2> xAndy_B) {
	for (auto i = 0u; i < 2; ++i) {
		for (auto j = 0u; j < xAndy_A.size(); ++j) {
			xAndy_B[i][j] += xAndy_A[i][j];
		}
	}
	return xAndy_B;
}

// goes from x to u using Zuhair's (4.21); returned vector has a list of all u+ 
// in order followed by a list of all u- in order
std::vector<char> UFromX(const std::string& x) {
	if (x.size() < 2) {
		std::cerr << "Error: asked to do exponent transform from X to U but "
			<< "there were only " << x.size() << " entries in X." << std::endl;
		return {};
	}

	std::vector<char> u(2*x.size() - 2);
	// the terms from x_1 through x_{n-1} are regular
	for (auto i = 0u; i < x.size()-1; ++i) {
		u[i] = 2*x[i];
		for (auto j = 0u; j < i; ++j) {
			u[x.size()-1 + j] += 2*x[i];
		}
	}
	// the last term, from x_n, is different
	for (auto j = 0u; j < x.size()-1; ++j) {
		u[x.size()-1 + j] += 2*x.back();
	}

	//std::cout << "Converted " << x << " into " << u << std::endl;
	return u;
}

namespace {
	std::unordered_map<std::string, std::vector<MatrixTerm_Final>> yCache;
} // anonymous namespace

std::vector<MatrixTerm_Final> ThetaFromY(const std::string y) {
	if (yCache.count(y) == 0) {
		std::vector<MatrixTerm_Intermediate> inter(YTildeFromY(y));
		std::vector<MatrixTerm_Final> ret(ThetaFromYTilde(inter));
		yCache.emplace(y, std::move(ret));
	}

	return yCache[y];
}

// convert from y to y-tilde following (4.26).
//
// this is by far the most intensive of the coordinate transformations: it has
// u biproducts in addition to the yTilde that you want, but worst of all it has
// two terms, so you end up with a sum of return terms, each with some 
// binomial-derived coefficient.
std::vector<MatrixTerm_Intermediate> YTildeFromY(const std::string& y) {
	//std::cout << "Transforming this y: " << y << std::endl;
	std::vector<MatrixTerm_Intermediate> ret;
	// i here is the i'th particle (pair); each one only sees those whose
	// numbers are lower, so we can treat them using lower particle numbers
	//
	// this loop goes to y.size()-1 because the last term, y_n, is special. The
	// first term is not special mathematically, but it's easier to code if we
	// do it separately as well

	// we begin with y_n because it's different from the others; it's restricted
	// to be the negative sum of all other y_i, so we can replace it directly at
	// the beginning; sadly this produces a bunch of terms
	std::vector<YTerm> yTerms = EliminateYn(y);
	//std::cout << "yTerms: " << std::endl;
	//for (const auto& yTerm : yTerms) std::cout << yTerm << std::endl;
	
	for (const YTerm& yTerm : yTerms) {
		std::vector<MatrixTerm_Intermediate> termsFromThisYTerm;
		// y_1 handled separately because it doesn't need multinomials. We 
		// always do this step even if yTerm[0] == 0 because it puts the
		// coefficient in and provides somewhere for the other terms to combine
		termsFromThisYTerm.emplace_back(1);
		termsFromThisYTerm.back().coefficient = yTerm.coeff;
		termsFromThisYTerm.back().uPlus[0] = yTerm[0];
		termsFromThisYTerm.back().uMinus[0] = yTerm[0];
		termsFromThisYTerm.back().yTilde[0] = yTerm[0];

		// terms between y_2 and y_{n-1}, inclusive (beware 1- vs 0-indexing)
		for (auto i = 1u; i < yTerm.size(); ++i) {
			if (yTerm[i] == 0) continue;
			std::vector<MatrixTerm_Intermediate> termsFromThisY;
			for (char l = 0; l <= yTerm[i]; ++l) {
				for (const auto& nAndm : Multinomial::GetMVectors(i, yTerm[i]-l)) {
					std::vector<MatrixTerm_Intermediate> newTerms(
							YTildeTerms(i, yTerm[i], l, nAndm) );
					termsFromThisY.insert(termsFromThisY.end(),
							newTerms.begin(), newTerms.end());
					//ret.insert(ret.end(), newTerms.begin(), newTerms.end());
				}
			}
			termsFromThisYTerm = MultiplyIntermediateTerms(
					termsFromThisYTerm, termsFromThisY );
		}
		//std::cout << "Transformed " << yTerm << " into these:" << std::endl;
		//for (auto& term : termsFromThisYTerm) std::cout << term << std::endl;
		ret.insert(ret.end(), termsFromThisYTerm.begin(),
				termsFromThisYTerm.end() );
	}

	// handle y_n separately here; the answer is actually wrong in Zuhair's
	// notes, so the below code is wrong as well
	/*Multinomial::Initialize(y.size()-1, y.back());
	for (char l = 0; l < y.back(); ++l) {
		for (const auto& nAndm : Multinomial::GetMVectors(y.size()-1, y.back() - l)) {
			MatrixTerm_Intermediate lastTerm(YTildeLastTerm(y.size() - 1, 
						y.back(), l, 
						std::vector<char>(nAndm.begin()+1, nAndm.end()) ) );
			ret.push_back(std::move(lastTerm));
		}
	}*/

	//std::cout << "Got these:" << std::endl;
	//for (auto& term : ret) std::cout << term << std::endl;

	return ret;
}

std::vector<YTerm> EliminateYn(const std::string& y) {
	std::vector<YTerm> output;
	for (auto nAndm : Multinomial::GetMVectors(y.size()-1, y.back())) {
		coeff_class coeff = Multinomial::Lookup(y.size()-1, nAndm);
		if (y.back() % 2 == 1) coeff = -coeff;
		do {
			output.emplace_back(coeff, y, nAndm);
		} while (std::prev_permutation(nAndm.begin()+1, nAndm.end()));
	}
	return output;
}

// part of the transformation in coordinates.pdf; mVectors must have length i,
// which accounts for the 0-indexing in C++ (it has length i-1 in the PDF).
/*MatrixTerm_Intermediate YTildeTerm(const unsigned int i, const char a, 
		const char l, const std::vector<char>& mVector) {
	MatrixTerm_Intermediate ret(i+1);

	for (auto j = 0u; j <= i-1; ++j) {
		ret.uPlus[j] = mVector[j];
		ret.yTilde[j] = mVector[j];
		ret.uMinus[j] = a;
		for (auto k = j+1; k < i; ++k) {
			ret.uMinus[j] += mVector[k];
		}
	}
	ret.uPlus[i] = 2*a - l;
	ret.uMinus[i] = l;
	ret.yTilde[i] = l;

	ret.coefficient = YTildeCoefficient(a, l, mVector);
	return ret;
}*/

// i is the y from y_i; a is the exponent of y_i; l is a binomial index; nAndm
// is a multinomial vector with total order a-l
std::vector<MatrixTerm_Intermediate> YTildeTerms(
		const unsigned int i, const char a, const char l, std::string nAndm) {
	std::vector<MatrixTerm_Intermediate> ret;

	coeff_class coeff = YTildeCoefficient(a, l, nAndm);
	do {
		ret.emplace_back(i+1);
		for (auto j = 0u; j <= i-1; ++j) {
			ret.back().uPlus[j] = nAndm[j+1];
			ret.back().yTilde[j] = nAndm[j+1];
			ret.back().uMinus[j] = a;
			// for (auto k = j+1; k < i; ++k) {
			for (auto k = 0u; k < j; ++k) {
				ret.back().uMinus[j] += nAndm[k+1];
				// ret.back().uMinus[k] += nAndm[j+1];
			}
		}
		ret.back().uPlus[i] = 2*a - l;
		ret.back().uMinus[i] = l;
		ret.back().yTilde[i] = l;

		ret.back().coefficient = coeff;
	} while (std::prev_permutation(nAndm.begin()+1, nAndm.end()));
	return ret;
}

std::vector<MatrixTerm_Intermediate> MultiplyIntermediateTerms(
		const std::vector<MatrixTerm_Intermediate>& termsA, 
		const std::vector<MatrixTerm_Intermediate>& termsB) {
	if (termsA.size() == 0 || termsB.size() == 0) {
		std::cerr << "Warning: asked to multiply two term lists but one was "
			<< "empty. A: " << termsA.size() << "; B: " << termsB.size() << "." 
			<< std::endl;
		return termsA.size() == 0 ? termsB : termsA;
	}
	//if (termsA.size() == 0) return termsB;
	//if (termsB.size() == 0) return termsA;
	std::vector<MatrixTerm_Intermediate> output;
	for (const auto& termA : termsA) {
		for (const auto& termB : termsB) {
			output.push_back(termA * termB);
		}
	}
	return output;
}

/*
// this is the same idea as the above, just for the special term with i=n
MatrixTerm_Intermediate YTildeLastTerm(const unsigned int n, const char a, 
		const char l, const std::vector<char>& mVector) {
	MatrixTerm_Intermediate ret(n-1);
	ret.uPlus = mVector;
	ret.yTilde = mVector;

	for (auto j = 1u; j <= n-1; ++j) {
		ret.uMinus[j-1] = a;
		for (auto k = j+1; k <= n-1; ++k) {
			ret.uMinus[j-1] += mVector[k-1];
		}
	}
	ret.uPlus[n-2] += l;
	ret.yTilde[n-2] += l;

	ret.coefficient = YTildeLastCoefficient(a, l, mVector);
	return ret;
}
*/

// the coefficient of a YTildeTerm, i.e. everything that's not a u or yTilde
coeff_class YTildeCoefficient(const char a, const char l, 
		const std::string& nAndm) {
	//coeff_class ret = ExactBinomial(a, l);
	coeff_class ret = Multinomial::Choose(2, a, {static_cast<char>(a-l),l});
	ret *= Multinomial::Lookup(nAndm.size()-1, nAndm);
	if ((a-l) % 2 == 1) ret = -ret;
	return ret;
}

/*
// The coefficient of a YTildeLastTerm, i.e. everything that's not a u or yTilde
// Basically the same as the YTildeCoefficient; only the sign is different
coeff_class YTildeLastCoefficient(const char a, const char l, 
		const std::vector<char>& mVector) {
	coeff_class ret = ExactBinomial(a, l);
	ret *= Multinomial::Choose(mVector.size(), a-l, mVector);
	if (a % 2 == 1) ret = -ret;
	return ret;
}
*/

// convert from y-tilde to sines and cosines of theta following (4.32).
//
// returned vector has sines of all components in order followed by all cosines
std::vector<MatrixTerm_Final> ThetaFromYTilde(
		std::vector<MatrixTerm_Intermediate>& intermediateTerms) {
	std::vector<MatrixTerm_Final> ret;

	for (auto& term : intermediateTerms) {
		// sine[i] appears in all yTilde[j] with j > i (strictly greater)
		std::vector<char> sines(term.yTilde.size()-1, 0);
		for (auto i = 0u; i < sines.size(); ++i) {
			for (auto j = i+1; j < term.yTilde.size(); ++j) {
				sines[i] += term.yTilde[j];
			}
		}
		ret.emplace_back(term.coefficient,
				std::move(term.uPlus), 
				std::move(term.uMinus),
				std::move(sines),
				std::move(term.yTilde) ); // this brings one spurious component
	}

	return ret;
}

std::vector<char> AddVectors(const std::vector<char>& A, 
		const std::vector<char>& B) {
	std::vector<char> output(std::max(A.size(), B.size()), 0);
	for (auto i = 0u; i < std::min(A.size(), B.size()); ++i) {
		output[i] = A[i] + B[i];
	}
	if (A.size() > B.size()) {
		for (auto i = B.size(); i < A.size(); ++i) {
			output[i] = A[i];
		}
	}
	if (B.size() > A.size()) {
		for (auto i = A.size(); i < B.size(); ++i) {
			output[i] = B[i];
		}
	}
	//std::cout << A << " + " << B << " = " << output << std::endl;
	return output;
	/* Better version:
	 * const std::vector<char>* aP = &A;
	 * const std::vector<char>* bP = &B;
	 * if (A.size() < B.size()) std::swap(aP, bP);
	 * std::vector<char> output(*aP);
	 * for (std::size_t i = 0; i < bP->size(); ++i) output[i] += (*bP)[i];
	 * return output;
	 */
}

/*std::vector<char> AddVectors(const std::vector<char>& A, 
		std::vector<char> B, const std::vector<std::size_t>& permutation) {
	std::vector<char> output(std::max(A.size(), B.size()), 0);
	// if this is slow, we should pad both A and B to full size before calling
	if (B.size() < permutation.size()) B.resize(permutation.size(), 0);

	for (auto i = 0u; i < std::min(A.size(), B.size()); ++i) {
		output[i] = A[i] + B[permutation[i]];
	}
	if (A.size() > B.size()) {
		for (auto i = B.size(); i < A.size(); ++i) {
			output[i] = A[i];
		}
	}
	if (B.size() > A.size()) {
		for (auto i = A.size(); i < B.size(); ++i) {
			output[i] = B[permutation[i]];
		}
	}
	std::cout << A << " + " << B << " = " << output << std::endl;
	return output;
}*/

// combines two u-and-theta coordinate wavefunctions (called F in Zuhair's
// notes), each corresponding to one monomial. Each is a sum over many terms, so
// combining them involves multiplying out two sums.
std::vector<MatrixTerm_Final> CombineTwoFs(const std::vector<MatrixTerm_Final>& F1,
		const std::vector<MatrixTerm_Final>& F2) {
	/*std::cout << "Combining two Fs with sizes " << F1.size() << " and " 
		<< F2.size() << std::endl;*/
	std::vector<MatrixTerm_Final> ret;
	for (auto& term1 : F1) {
		for (auto& term2 : F2) {
			ret.emplace_back(term1.coefficient * term2.coefficient,
					AddVectors(term1.uPlus, term2.uPlus),
					AddVectors(term1.uMinus, term2.uMinus),
					AddVectors(term1.sinTheta, term2.sinTheta),
					AddVectors(term1.cosTheta, term2.cosTheta));
		}
	}
	//std::cout << "Combined Fs has " << ret.size() << " terms." << std::endl;
	return ret;
}

// WARNING: if type == MAT_MASS this breaks the MatrixTerm_Final vector fed into
// it by changing each term's uMinus entries by -2
coeff_class FinalResult(std::vector<MatrixTerm_Final>& exponents,
		const MATRIX_TYPE type) {
	if (exponents.size() == 0) {
		std::cout << "No exponents detected; returning 1." << std::endl;
		return 1;
	}
	// auto n = exponents.front().uPlus.size() + 1;
	coeff_class totalFromIntegrals = 0;
	for (auto& term : exponents) {
		if (type == MAT_INNER) {
			// just do the integrals
			totalFromIntegrals += DoAllIntegrals(term);
		} else if (type == MAT_MASS) {
			// sum over integral results for every possible 1/x
			term.uPlus[0] -= 2;
			totalFromIntegrals += DoAllIntegrals(term);
			for (std::size_t i = 1; i < term.uPlus.size(); ++i) {
				term.uPlus[i-1] += 2;
				term.uMinus[i-1] -= 2;
				term.uPlus[i] -= 2;
				totalFromIntegrals += DoAllIntegrals(term);
			}
			term.uPlus.back() += 2;
			term.uMinus.back() -= 2;
			totalFromIntegrals += DoAllIntegrals(term);
		}
	}
	/*std::cout << "Returning FinalResult = " << InnerProductPrefactor(n) << " * " 
		<< totalFromIntegrals << std::endl;
	return InnerProductPrefactor(n) * totalFromIntegrals;*/
	//std::cout << "Returning FinalResult = " << MassMatrixPrefactor(n) << " * " 
		//<< totalFromIntegrals << std::endl;
	// if (type == MAT_INNER) {
		return totalFromIntegrals;
	// } else if (type == MAT_MASS) {
		// return MassMatrixPrefactor(n) * totalFromIntegrals;
	// }
}

coeff_class InnerProductPrefactor(const char n) {
	return MassMatrixPrefactor(n)/n;
}

coeff_class MassMatrixPrefactor(const char n) {
	coeff_class denominator = std::tgamma(n); // tgamma is the "true" gamma fcn
	denominator *= std::pow(16, n-1);
	denominator *= std::pow(M_PI, 2*n-3);
	//std::cout << "PREFACTOR: " << 2/denominator << std::endl;
	return 2/denominator;
}

// integrals ------------------------------------------------------------------

// cache results using hash tables keyed by a and b
namespace {
	std::unordered_map<std::array<coeff_class,2>, coeff_class,
		boost::hash<std::array<coeff_class,2>> > uCache;
	std::unordered_map<std::array<coeff_class,2>, coeff_class,
		boost::hash<std::array<coeff_class,2>> > thetaCache;
} // anonymous namespace

coeff_class DoAllIntegrals(const MatrixTerm_Final& term) {
	std::size_t n = term.uPlus.size() + 1;
	coeff_class output = term.coefficient;

	// do the u integrals first
	for (auto i = 0u; i < n-1; ++i) {
		output *= UIntegral(term.uPlus[i] + 3, 5*(n - i) - 7 + term.uMinus[i]);
	}

	// now the theta integrals; sineTheta.size() is n-2, but cosTheta.size()
	// is n-1 with the last component being meaningless. All but the last 
	// one are short, while the last one is long
	//
	// these have constant terms which differ from Zuhair's because his i
	// starts at 1 instead of 0
	if (n >= 3) {
		for (auto i = 0u; i < n-3; ++i) {
			output *= ThetaIntegral_Short(n-3 - i + term.sinTheta[i],
					term.cosTheta[i] );
		}
		output *= ThetaIntegral_Long(term.sinTheta[n-3], term.cosTheta[n-3]);
	}
	// std::cout << term.coefficient << " * {" << term.uPlus << ", " << term.uMinus
		// << ", " << term.sinTheta << ", " << term.cosTheta << "} -> " << output 
		// << std::endl;
	return output;
}

// this is the integral over the "u" variables; it uses a hypergeometric 
// identity to turn 2F1(a, b, c, -1) -> 2^(-a)*2F1(a, c-b, c, 1/2)
//
// follows the conventions of Zuhair's 5.34; a is the exponent of u_i+ and
// b is the exponent of u_i-
coeff_class UIntegral(const coeff_class a, const coeff_class b) {
	//std::cout << "UIntegral(" << a << ", " << b << ")" << std::endl;
	std::array<coeff_class,2> abArray{{a,b}};
	if (b < a) std::swap(abArray[0], abArray[1]);
	if (uCache.count(abArray) == 1) return uCache.at(abArray);

	coeff_class ret = gsl_sf_hyperg_2F1(1, (a+b)/2 + 2, b/2 + 2, 0.5)/(b + 2);
	ret += gsl_sf_hyperg_2F1(1, (a+b)/2 + 2, a/2 + 2, 0.5)/(a + 2);
	ret *= std::pow(std::sqrt(2), -(a+b));
	uCache.emplace(abArray, ret);
	return ret;
}

// this is the integral over the "theta" veriables from 0 to pi; it implements 
// Zuhair's 5.35, where a is the exponent of sin(theta) and b is the exponent of 
// cos(theta).
//
// results are cached by (a,b); since a and b are symmetric, we only store the
// results with a <= b, swapping the two parameters if they're the other order
coeff_class ThetaIntegral_Short(const coeff_class a, const coeff_class b) {
	if (static_cast<int>(b) % 2 == 1) return 0;
	std::array<coeff_class,2> abArray{{a,b}};
	if (b < a) std::swap(abArray[0], abArray[1]);
	if (thetaCache.count(abArray) == 1) return thetaCache.at(abArray);

	coeff_class ret = std::exp(std::lgamma((1+a)/2) + std::lgamma((1+b)/2) 
			- std::lgamma((2 + a + b)/2) );
	thetaCache.emplace(abArray, ret);
	return ret;
}

// this is the integral over the "theta" veriables from 0 to 2pi; it implements 
// Zuhair's 5.36, where a is the exponent of sin(theta) and b is the exponent of
// cos(theta).
coeff_class ThetaIntegral_Long(const coeff_class a, const coeff_class b) {
	if (static_cast<int>(a) % 2 == 1) return 0;
	return 2*ThetaIntegral_Short(a, b);
}

} // namespace MatrixInternal
