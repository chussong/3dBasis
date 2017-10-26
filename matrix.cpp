#include "matrix.hpp"

DMatrix MassMatrix(const Basis<Mono>& basis) {
	DMatrix massMatrix(basis.size(), basis.size());
	for (std::size_t i = 0; i < basis.size(); ++i) {
		massMatrix(i, i) = MatrixInternal::MassMatrixTerm(basis[i], basis[i]);
		for (std::size_t j = i+1; j < basis.size(); ++j) {
			massMatrix(i, j) = MatrixInternal::MassMatrixTerm(basis[i], basis[j]);
			massMatrix(j, i) = massMatrix(i, j);
		}
	}
	return massMatrix;
}

namespace MatrixInternal {

MatrixTerm_Intermediate::MatrixTerm_Intermediate(const size_t n): coefficient(1),
	uPlus(n), uMinus(n), yTilde(n) {
}

void MatrixTerm_Intermediate::Resize(const size_t n) {
	uPlus.resize(n);
	uMinus.resize(n);
	yTilde.resize(n);
}

// this n would be called (n-1) in Zuhair's equations
MatrixTerm_Final::MatrixTerm_Final(const size_t n): coefficient(1), uPlus(n), 
	uMinus(n), sinTheta(n-1), cosTheta(n-1) {
}

// maybe should be rvalue references instead? hopefully it's the same
MatrixTerm_Final::MatrixTerm_Final(const coeff_class coefficient, const std::vector<char>& uPlus, 
		const std::vector<char>& uMinus, const std::vector<char>& sinTheta,
		const std::vector<char>& cosTheta): coefficient(coefficient), 
	uPlus(uPlus), uMinus(uMinus), sinTheta(sinTheta), cosTheta(cosTheta) {
}

// notice that there are fewer Thetas than Us, so the vectors aren't all the
// same size; also, the n here would be called (n-1) in Zuhair's equations
void MatrixTerm_Final::Resize(const size_t n) {
	uPlus.resize(n);
	uMinus.resize(n);
	sinTheta.resize(n-1);
	cosTheta.resize(n-1);
}

coeff_class MassMatrixTerm(const Mono& A, const Mono& B) {
	std::vector<MatrixTerm_Final> fFromA(MatrixTermsFromMono(A));
	std::vector<MatrixTerm_Final> fFromB(MatrixTermsFromMono(B));
	// the below step clears the two original vectors automatically
	std::vector<MatrixTerm_Final> combinedFs(CombineTwoFs(fFromA, fFromB));

	return FinalResult(combinedFs);
}

// takes a Mono and gets the exponents of the u and theta variables * some coeff
std::vector<MatrixTerm_Final> MatrixTermsFromMono(const Mono& input) {
	std::array<std::vector<char>,2> xAndy(ExponentExtractXY(input));
	std::vector<char> uFromX(ExponentUFromX(xAndy[0]));
	std::vector<MatrixTerm_Intermediate> inter(ExponentYTildeFromY(xAndy[1]));
	std::vector<MatrixTerm_Final> ret(ExponentThetaFromYTilde(inter));
	// u has contributions from both x and y, so we have to combine them
	for (auto& term : ret) {
		for (auto i = 0u; i < term.uPlus.size(); ++i) {
			term.uPlus[i] += uFromX[i];
			term.uMinus[i] += uFromX[term.uPlus.size() + i];
		}
	}
	return ret;
}

/*
// variable changes / coordinate transformations ------------------------------
//
// These are ways of computing a given value in (u, theta) given values in (x,y);
// I wrote them to get a feel for the coorinate system, and they're not really useful.
//
// These transformations generally have terms with coefficients that are partial
// products of lower-numbered elements in various vectors. The algorithms below
// take advantage of this by computing them as running products which are 
// updated progressively as the computation proceeds, meaning old elements 
// generally do not ever need to be re-accessed.

// extracts the x and y from a monomial (conventions from Zuhair's (4.15))
std::array<std::vector<coeff_class>, 2> ExtractXY(const Mono& extractFromThis) {
	coeff_class totalPm = extractFromThis.TotalPm();
	coeff_class mu = 1; // should be something else, but whatever

	// it'd be faster to preallocate the array with the right size and fill it
	// in, but I can't imagine this will actually take very long anyway
	std::vector<coeff_class> x, y;
	for (const auto& part : extractFromThis) {
		x.push_back(coeff_class{part.pm} / totalPm);
		y.push_back(coeff_class{part.pt} / mu);
	}
	return {x, y};
}

// transforms from x variables to u variables according to Zuhair's (4.21, 4.22)
std::vector<coeff_class> UFromX(const std::vector<coeff_class>& x) {
	if (x.size() < 1) {
		std::cerr << "Error: asked to convert an empty x vector to u variables."
			<< std::endl;
		return {};
	}
	std::vector<coeff_class> u(x.size() - 1);

	coeff_class coefficientOfX = 2;
	for (auto i = 0u; i < u.size(); ++i) {
		u[i] = (coefficientOfX * x[i]) - 1;
		coefficientOfX *= 2/(1-u[i]);
	}
	return u;
}

// transforms from y variables to y-tilde variables following (4.26, 4.27)
std::vector<coeff_class> YTildeFromYAndU(const std::vector<coeff_class>& y,
		const std::vector<coeff_class>& u) {
	std::vector<coeff_class> yTilde(u.size());
	coeff_class productOfMinuses = 1;
	std::vector<coeff_class> secondTermPieces;
	for (auto i = 0u; i < yTilde.size(); ++i) {
		coeff_class uMinus = std::sqrt((1-u[i])/2);
		coeff_class uPlus = std::sqrt((1+u[i])/2);

		yTilde[i] = y[i]/(productOfMinuses * uMinus*uPlus);
		coeff_class secondTerm = 0;
		for (auto& piece : secondTermPieces) {
			secondTerm += piece;
			piece *= uMinus; // this prepares the piece for the next i iteration
		}
		yTilde[i] += secondTerm * uPlus / uMinus;

		// more preparation for the next iteration
		secondTermPieces.push_back(yTilde[i] * uPlus);
		productOfMinuses *= uMinus;
	}
	return yTilde;
}

// transforms from y-tilde variables to cosines of theta following (4.32)
std::vector<coeff_class> SinThetaFromYTilde(const std::vector<coeff_class>& yT){
	coeff_class productOfSines = 1;
	std::vector<coeff_class> sinTheta(yT.size() - 1);

	for (auto i = 0u; i < sinTheta.size(); ++i) {
		coeff_class cosTheta = yT[i]/productOfSines;
		sinTheta[i] = std::sqrt(1 - cosTheta*cosTheta);
		productOfSines *= sinTheta[i];
	}
	return sinTheta;
}*/

// exponent transformations ---------------------------------------------------
//
// rather than transforming the variables themselves, these take in a list of 
// exponents in one coordinate system and give out lists of exponents in the new 
// coordinate system.

// just gets the exponents of each x and y, so we don't have to worry about
// things like dividing by P or \mu
std::array<std::vector<char>,2> ExponentExtractXY(const Mono& extractFromThis) {
	std::vector<char> x, y;
	for (auto i = 0u; i < extractFromThis.NParticles(); ++i) {
		x.push_back(extractFromThis.Pm(i));
		y.push_back(extractFromThis.Pm(i));
	}
	return {{x, y}};
}

// goes from x to u using Zuhair's (4.21); returned vector has a list of all u+ 
// in order followed by a list of all u- in order
std::vector<char> ExponentUFromX(const std::vector<char>& x) {
	if (x.size() < 2) {
		std::cerr << "Error: asked to do exponent transform from X to U but "
			<< "there were only " << x.size() << " entries in X." << std::endl;
		return {};
	}

	std::vector<char> u(2*x.size() - 2);
	for (auto i = 0u; i < x.size()-1; ++i) {
		u[i] = 2*x[i];
		for (auto j = 0u; j < i; ++j) {
			u[x.size()-1 + j] += 2*x[i];
		}
	}
	for (auto j = 0u; j < x.size()-1; ++j) {
		u[x.size()-1 + j] += 2*x.back();
	}

	return u;
}

// convert from y to y-tilde following (4.26).
//
// this is by far the most intensive of the coordinate transformations: it has
// u biproducts in addition to the yTilde that you want, but worst of all it has
// two terms, so you end up with a sum of return terms, each with some 
// binomial-derived coefficient.
std::vector<MatrixTerm_Intermediate> ExponentYTildeFromY(const std::vector<char>& y) {
	std::vector<MatrixTerm_Intermediate> ret;
	// i here is the i'th particle (pair); each one only sees those whose
	// numbers are lower, so we can treat them using lower particle numbers
	//
	// this loop goes to y.size()-1 because the last term, y_n, is special
	for (auto i = 0u; i < y.size()-1; ++i) {
		Multinomial::Initialize(i, y[i]);
		for (char l = 0; l <= y[i]; ++l) {
			for (const auto& nAndm : Multinomial::GetMVectors(y[i]-l)) {
				MatrixTerm_Intermediate newTerm(YTildeTerm(i, y[i], l, 
							std::vector<char>(nAndm.begin()+1, nAndm.end())));
				ret.push_back(std::move(newTerm));
			}
		}
	}
	// handle y_n separately here
	Multinomial::Initialize(y.size()-1, y.back());
	for (char l = 0; l < y.back(); ++l) {
		for (const auto& nAndm : Multinomial::GetMVectors(y.back() - l)) {
			MatrixTerm_Intermediate lastTerm(YTildeLastTerm(y.size() - 1, 
						y.back(), l, 
						std::vector<char>(nAndm.begin()+1, nAndm.end()) ) );
			ret.push_back(std::move(lastTerm));
		}
	}
	return ret;
}

// part of the transformation in coordinates.pdf; mVectors must have length i,
// which accounts for the 0-indexing in C++ (it has length i-1 in the PDF).
MatrixTerm_Intermediate YTildeTerm(const unsigned int i, const char a, 
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
	ret.uPlus[i] = a;
	ret.uMinus[i] = l;
	ret.yTilde[i] = l;

	ret.coefficient = YTildeCoefficient(a, l, mVector);
	return ret;
}

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

// the coefficient of a YTildeTerm, i.e. everything that's not a u or yTilde
coeff_class YTildeCoefficient(const char a, const char l, 
		const std::vector<char>& mVector) {
	coeff_class ret = ExactBinomial(a, l);
	ret *= Multinomial::Choose(a-l, mVector);
	if (a-l % 2 == 1) ret = -ret;
	return ret;
}

// The coefficient of a YTildeLastTerm, i.e. everything that's not a u or yTilde
// Basically the same as the YTildeCoefficient; only the sign is different
coeff_class YTildeLastCoefficient(const char a, const char l, 
		const std::vector<char>& mVector) {
	coeff_class ret = ExactBinomial(a, l);
	ret *= Multinomial::Choose(a-l, mVector);
	if (a % 2 == 1) ret = -ret;
	return ret;
}

// convert from y-tilde to sines and cosines of theta following (4.32).
//
// returned vector has sines of all components in order followed by all cosines
std::vector<MatrixTerm_Final> ExponentThetaFromYTilde(
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
	return output;
}

// warning: when combining F1 and F2, this clears the original ones
std::vector<MatrixTerm_Final> CombineTwoFs(std::vector<MatrixTerm_Final>& F1,
		std::vector<MatrixTerm_Final>& F2) {
	std::vector<MatrixTerm_Final> ret;
	for (auto& term1 : F1) {
		for (auto& term2 : F2) {
			// the below probably won't compile, and we'll have to define 
			// a specialized Plus(std::vector<char>, std::vector<char>)
			ret.emplace_back(term1.coefficient * term2.coefficient,
					AddVectors(term1.uPlus, term2.uPlus),
					AddVectors(term1.uMinus, term2.uMinus),
					AddVectors(term1.sinTheta, term2.sinTheta),
					AddVectors(term1.cosTheta, term2.cosTheta) );
		}
	}
	F1.clear();
	F2.clear();
	return ret;
}

coeff_class FinalResult(const std::vector<MatrixTerm_Final>& exponents) {
	if (exponents.size() == 0) return MatrixPrefactor(2);
	auto n = exponents.front().uPlus.size() + 1;
	coeff_class integral = 1;
	for (const auto& term : exponents) {
		// do the u integrals first; uPlus.size() is n-1
		for (auto i = 0u; i < n-1; ++i) {
			integral *= UIntegral(3 + term.uPlus[i], 
					5*(n - i) - 2 + term.uMinus[i]);
		}

		// now the theta integrals; sineTheta.size() is n-2, but cosTheta.size()
		// is n-1 with the last component being meaningless
		for (auto i = 0u; i < n-3; ++i) {
			integral *= ThetaIntegral_Short(n - 2 - i + term.sinTheta[i],
					term.cosTheta[i] );
		}
		integral *= ThetaIntegral_Long(term.sinTheta[n-3], term.cosTheta[n-3]);
	}
	return MatrixPrefactor(n) * integral;
}

coeff_class MatrixPrefactor(const char n) {
	coeff_class denominator = std::tgamma(n); // tgamma is the "true" gamma fcn
	denominator *= std::pow(16, n-1);
	denominator *= std::pow(M_PI, 2*n-3);
	// is m^2 some kind of mass? // ret *= m*m;
	return 1/denominator;
}

// integrals ------------------------------------------------------------------

// cache results using hash tables keyed by a and b
namespace {
	std::unordered_map<std::array<coeff_class,2>, coeff_class,
		boost::hash<std::array<coeff_class,2>> > uCache;
	std::unordered_map<std::array<coeff_class,2>, coeff_class,
		boost::hash<std::array<coeff_class,2>> > thetaCache;
} // anonymous namespace

// this is the integral over the "u" variables; it uses a hypergeometric 
// identity to turn 2F1(a, b, c, -1) -> 2^(-a)*2F1(a, c-b, c, 1/2)
//
// follows the conventions of Zuhair's 4.34; a is the exponent of u_i+ and
// b is the exponent of u_i-
coeff_class UIntegral(const coeff_class a, const coeff_class b) {
	std::array<coeff_class,2> abArray;
	if (uCache.count(abArray) == 1) return uCache.at(abArray);

	coeff_class ret = gsl_sf_hyperg_2F1(1, (a+b)/2 + 2, b/2 + 2, 0.5)/(b + 2);
	ret += gsl_sf_hyperg_2F1(1, (a+b)/2 + 2, a/2 + 2, 0.5)/(a + 2);
	ret *= std::pow(std::sqrt(2), -(a+b));
	uCache.emplace(abArray, ret);
	return ret;
}

// this is the integral over the "theta" veriables from 0 to pi; it implements 
// Zuhair's 4.35, where a is the exponent of sin(theta) and b is the exponent of 
// cos(theta).
coeff_class ThetaIntegral_Short(const coeff_class a, const coeff_class b) {
	if (static_cast<int>(b) % 2 == 1) return 0;
	std::array<coeff_class,2> abArray;
	if (thetaCache.count(abArray) == 1) return thetaCache[abArray];

	coeff_class ret = std::exp(std::lgamma((1+a)/2) + std::lgamma((1+b)/2) 
			- std::lgamma((2 + a + b)/2) );
	thetaCache.emplace(abArray, ret);
	return ret;
}

// this is the integral over the "theta" veriables from 0 to 2pi; it implements 
// Zuhair's 4.36, where a is the exponent of sin(theta) and b is the exponent of
// cos(theta).
coeff_class ThetaIntegral_Long(const coeff_class a, const coeff_class b) {
	if (static_cast<int>(a) % 2 == 1 || static_cast<int>(b) % 2 == 1) return 0;
	return 2*ThetaIntegral_Short(a, b);
}

} // namespace Matrix
