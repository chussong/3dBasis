#include "matrix.hpp"

// Fock space part (ONLY) of the inner product between two monomials
coeff_class InnerFock(const Mono& A, const Mono& B) {
    //B.ChangePm(0, B.Pm(0)+1); // this will be reverted in MatrixTerm below
    return MatrixInternal::MatrixTerm(A, B, MAT_INNER);
    //return MatrixInternal::MatrixTerm(A, B, MatrixInter::MAT_INNER)/A.NParticles();
}

// inner product between two partitions of monomials
coeff_class InnerProduct(const Mono& A, const Mono& B) {
    // DVector muIntegrals = MuIntegral(A, B, partitions, partWidth, MAT_INNER);
    return MatrixInternal::MatrixTerm(A, B, MAT_INNER);
        // *muIntegrals(part);
}

// creates a gram matrix for the given basis using the Fock space inner product
//
// this returns the rank 2 matrix containing only the Fock part of the product
DMatrix GramFock(const Basis<Mono>& basis) {
    return MatrixInternal::Matrix(basis, 0, 0, MAT_INNER);
}

// creates a gram matrix for the given basis using the Fock space inner product
// 
// this returns the rank 4 tensor relating states with different partitions
DMatrix GramMatrix(const Basis<Mono>& basis, const std::size_t partitions, 
        const coeff_class partWidth) {
    return MatrixInternal::Matrix(basis, partitions, partWidth, MAT_INNER);
}

// creates a mass matrix M for the given monomials. To get the mass matrix of a 
// basis of primary operators, one must express the primaries as a matrix of 
// vectors, A, and multiply A^T M A.
DMatrix MassMatrix(const Basis<Mono>& basis, const std::size_t partitions, 
        const coeff_class partWidth) {
    return MatrixInternal::Matrix(basis, partitions, partWidth, MAT_MASS);
}

DMatrix KineticMatrix(const Basis<Mono>& basis, const std::size_t partitions,
        const coeff_class partWidth) {
    return MatrixInternal::Matrix(basis, partitions, partWidth, MAT_KINETIC);
}

// creates a matrix of n->n interactions between the given basis's monomials
DMatrix InteractionMatrix(const Basis<Mono>& basis, const std::size_t partitions, 
        const coeff_class partWidth) {
    return MatrixInternal::Matrix(basis, partitions, partWidth, MAT_INTER_SAME_N);
}

namespace MatrixInternal {

// static hash tables for memoizing slow steps
namespace {
    // all matrices: map from {x,y}->{u,yTilde}
    std::unordered_map<std::string, std::vector<MatrixTerm_Intermediate>>
        intermediateCache;
    // direct matrices: map from {x,y}->{u,theta}
    std::unordered_map<std::string, std::vector<MatrixTerm_Final>> directCache;
    // interaction matrices: map from {x,y}->{u,r,theta}
    std::unordered_map<std::string, std::vector<InteractionTerm_Step2>>
        interactionCache;
    // interaction n+2 matrices: map from {x,y}->{u,theta}
    std::unordered_map<std::string, std::vector<MatrixTerm_Final>> nPlus2Cache;

    // integrals: map from (a,b)->#
	std::unordered_map<std::array<builtin_class,2>, builtin_class,
		boost::hash<std::array<builtin_class,2>> > uCache;
	std::unordered_map<std::array<builtin_class,2>, builtin_class,
		boost::hash<std::array<builtin_class,2>> > thetaCache;
} // anonymous namespace

YTerm::YTerm(const coeff_class coeff, const std::string& y, 
		const std::string& nAndm): coeff(coeff), y(y.begin(), y.end()-1) {
	for (std::size_t i = 0; i < size(); ++i) {
		this->y[i] += nAndm[i+1];
	}
}

std::ostream& operator<<(std::ostream& os, const YTerm& out) {
	return os << out.coeff << " * " << out.y;
}

MatrixTerm_Intermediate::MatrixTerm_Intermediate(const size_t n): coeff(1),
	uPlus(n), uMinus(n), yTilde(n) {
}

void MatrixTerm_Intermediate::Resize(const size_t n) {
    uPlus.resize(n);
    uMinus.resize(n);
    yTilde.resize(n);
}

MatrixTerm_Intermediate operator*(
		const MatrixTerm_Intermediate& A, MatrixTerm_Intermediate B) {
	B.coeff *= A.coeff;
	B.uPlus  = AddVectors(A.uPlus , B.uPlus );
	B.uMinus = AddVectors(A.uMinus, B.uMinus);
	B.yTilde = AddVectors(A.yTilde, B.yTilde);
	return B;
}

std::ostream& operator<<(std::ostream& os, const MatrixTerm_Intermediate& out) {
	return os << out.coeff << " * {" << out.uPlus << ", "
		<< out.uMinus << ", " << out.yTilde << "}";
}

// this n would be called (n-1) in Zuhair's equations
MatrixTerm_Final::MatrixTerm_Final(const size_t n): coeff(1), uPlus(n), 
	uMinus(n), sinTheta(n-1), cosTheta(n-1) {
}

// maybe should be rvalue references instead? hopefully it's the same
MatrixTerm_Final::MatrixTerm_Final(const coeff_class coeff, 
		const std::vector<char>& uPlus, const std::vector<char>& uMinus, 
		const std::vector<char>& sinTheta, const std::vector<char>& cosTheta): 
	coeff(coeff), uPlus(uPlus), uMinus(uMinus), sinTheta(sinTheta), 
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

std::ostream& operator<<(std::ostream& os, const InteractionTerm_Step2& out) {
    return os << out.coeff << " * {" << out.u << ", "
        << out.theta << ", " << out.r << "}";
}

// generically return direct or interaction matrix of the specified type
//
// NOTE: we will have to do discretizations before this function returns; in
// fact, we will have to do them inside MatrixTerm
DMatrix Matrix(const Basis<Mono>& basis, const std::size_t partitions, 
        const coeff_class partWidth, const MATRIX_TYPE type) {
    // partitions == 0 means that the Fock part has been requested by itself
    if (partitions == 0) {
        DMatrix fockPart(basis.size(), basis.size());
        for (std::size_t i = 0; i < basis.size(); ++i) {
            fockPart(i, i) = MatrixTerm(basis[i], basis[i], type);
            for (std::size_t j = i+1; j < basis.size(); ++j) {
                fockPart(i, j) = MatrixTerm(basis[i], basis[j], type);
                fockPart(j, i) = fockPart(i, j);
            }
        }

        return fockPart;
    } else {
        DMatrix output(basis.size()*partitions, basis.size()*partitions);
        for (std::size_t i = 0; i < basis.size(); ++i) {
            for (std::size_t j = 0; j < basis.size(); ++j) {
                output.block(i*partitions, j*partitions, partitions, partitions) = 
                    MatrixBlock(basis[i], basis[j], type, partitions, partWidth);
            }
        }
        return output;
    }
}

coeff_class MatrixTerm(const Mono& A, const Mono& B, const MATRIX_TYPE type) {
    if (type == MAT_INNER || type == MAT_MASS) {
        return MatrixTerm_Direct(A, B, type);
    } else if (type == MAT_KINETIC) {
        return MatrixTerm_Direct(A, B, MAT_INNER);
    } else if (type == MAT_INTER_N_PLUS_2) {
        return MatrixTerm_NPlus2(A, B);
    } else if (type == MAT_INTER_SAME_N) {
        throw std::logic_error("MatrixTerm: n-n interaction can't give scalar");
    } else {
        throw std::logic_error("MatrixTerm: unrecognized matrix type");
    }
}

DMatrix MatrixBlock(const Mono& A, const Mono& B, const MATRIX_TYPE type,
        const std::size_t partitions, const coeff_class partWidth) {
    if (type == MAT_INTER_SAME_N) {
        auto terms = MatrixTerm_Inter(A, B);
        DMatrix output = DMatrix::Zero(partitions, partitions);
        for (const auto& term : terms) {
            // std::cout << "Term: " << term.r << std::endl;
            output += term.coeff*MuPart(term.r, partitions, partWidth);
        }
        return output;
    } else {
        return MatrixTerm(A, B, type)*MuPart(A, B, partitions, partWidth, type);
    }
}

coeff_class MatrixTerm_Direct(const Mono& A, const Mono& B, const MATRIX_TYPE type) {
    //std::cout << "TERM: " << A.HumanReadable() << " x " << B.HumanReadable() 
            //<< std::endl;

    // degeneracy factors result from turning the ordered monomials into 
    // symmetric polynomials
    coeff_class degeneracy = 1;
    degeneracy *= Factorial(A.NParticles());
    // for (auto& count : A.CountIdentical()) degeneracy *= Factorial(count);
    for (auto& count : B.CountIdentical()) degeneracy *= Factorial(count);

    coeff_class prefactor = degeneracy*A.Coeff()*B.Coeff()*Prefactor(A, B, type);

    std::string xAndy_A = ExtractXY(A);
    std::string xAndy_B = ExtractXY(B);
    std::vector<MatrixTerm_Final> fFromA, fFromB, combinedFs;

    coeff_class total = 0;
    // do {
    fFromA = DirectTermsFromXY(xAndy_A);
    do {
        fFromB = DirectTermsFromXY(xAndy_B);
        combinedFs = CombineTwoFs(fFromA, fFromB);
        total += FinalResult(combinedFs, type);
    } while (PermuteXY(xAndy_B));
    // } while (PermuteXY(xAndy_A));

    return prefactor*total;
}

// the "inter" here means interacting, not intermediate, which is a bad name!
std::vector<InteractionTerm_Output> MatrixTerm_Inter(const Mono& A, 
        const Mono& B) {
    //std::cout << "INTERACTION: " << A.HumanReadable() << " x " 
    //<< B.HumanReadable() << std::endl;

    // degeneracy factors result from turning the ordered monomials into 
    // symmetric polynomials
    coeff_class degeneracy = 1;
    // degeneracy *= Factorial(A.NParticles());
    // degeneracy *= Factorial(B.NParticles());
    for (auto& count : A.CountIdentical()) degeneracy *= Factorial(count);
    for (auto& count : B.CountIdentical()) degeneracy *= Factorial(count);

    coeff_class prefactor = degeneracy*A.Coeff()*B.Coeff()*
        Prefactor(A, B, MAT_INTER_SAME_N);

    std::string xAndy_A = ExtractXY(A);
    std::string xAndy_B = ExtractXY(B);
    std::vector<MatrixTerm_Intermediate> fFromA, fFromB;
	std::vector<InteractionTerm_Step2> combinedFs;
    std::vector<InteractionTerm_Output> output;
    do {
        fFromA = InteractionTermsFromXY(xAndy_A);
        do {
            fFromB = InteractionTermsFromXY(xAndy_B);
            combinedFs = CombineInteractionFs(fFromA, fFromB);
            auto newTerms = InteractionOutput(combinedFs, MAT_INTER_SAME_N, 
                    prefactor);
            output.insert(output.end(), newTerms.begin(), newTerms.end());
        } while (PermuteXY(xAndy_B));
    } while (PermuteXY(xAndy_A));
    
    return output;
}

coeff_class MatrixTerm_NPlus2(const Mono& A, const Mono& B) {
    //std::cout << "N+2 TERM: " << A << " x " << B << std::endl;

    if (B.NParticles() != A.NParticles() + 2) {
        throw std::logic_error("MatrixTerm_NPlus2: n_B != n_A + 2");
    }

    // degeneracy factors result from turning the ordered monomials into 
    // symmetric polynomials
    coeff_class degeneracy = 1;
    // degeneracy *= Factorial(A.NParticles());
    for (auto& count : A.CountIdentical()) degeneracy *= Factorial(count);
    for (auto& count : B.CountIdentical()) degeneracy *= Factorial(count);

    coeff_class prefactor = degeneracy*A.Coeff()*B.Coeff()
        *Prefactor(A, B, MAT_INTER_N_PLUS_2);

    std::string xAndy_A = ExtractXY(A);
    std::string xAndy_B = ExtractXY(B);
    std::vector<MatrixTerm_Final> fFromA, fFromB, combinedFs;

    coeff_class total = 0;
    do {
        fFromA = DirectTermsFromXY(xAndy_A);
        do {
            fFromB = TermsFromXY_NPlus2(xAndy_B);
            combinedFs = CombineTwoFs(fFromA, fFromB);
            total += FinalResult_NPlus2(combinedFs);
        } while (PermuteXY(xAndy_B));
    } while (PermuteXY(xAndy_A));

    return prefactor*total;
}

// custom std::next_permutation for xAndy using particle precedence
bool PermuteXY(std::string& xAndy) {
    if (xAndy.size() % 2 != 0) {
        throw std::logic_error("Odd-length vector passed to PermuteXY");
    }
    if (xAndy.size() <= 2) return false;

    auto half = xAndy.size()/2;
    auto i = half - 1;

    while (i > 0) {
        auto i1 = i;
        i -= 1;
        if ( (xAndy[i] > xAndy[i1])
            || (xAndy[i] == xAndy[i1] && xAndy[half + i] > xAndy[half + i1]) ) {
            auto i2 = half;
            do {
                i2 -= 1;
            } while ( (xAndy[i] < xAndy[i2]) || 
                    (xAndy[i] == xAndy[i2] && xAndy[half+i] <= xAndy[half+i2]));
            std::swap(xAndy[i], xAndy[i2]);
            std::swap(xAndy[i + half], xAndy[i2 + half]);
            std::reverse(xAndy.begin() + i1, xAndy.begin() + half);
            std::reverse(xAndy.begin() + half + i1, xAndy.end());
            return true;
        }
    }
    std::reverse(xAndy.begin(), xAndy.begin() + half);
    std::reverse(xAndy.begin() + half, xAndy.end());
    return false;
}

const std::vector<MatrixTerm_Final>& DirectTermsFromXY(const std::string& xAndy)
{
    if (directCache.count(xAndy) == 0) {
        // copy so we can break it in the next function
        std::vector<MatrixTerm_Intermediate> intermediate 
            = InteractionTermsFromXY(xAndy);
        directCache.emplace(xAndy, ThetaFromYTilde(intermediate));
    }

    return directCache[xAndy];
}

const std::vector<MatrixTerm_Final>& TermsFromXY_NPlus2(
        const std::string& xAndy) {
    if (nPlus2Cache.count(xAndy) == 0) {
        std::vector<MatrixTerm_Intermediate> intermediate
            = InteractionTermsFromXY(xAndy);
        nPlus2Cache.emplace(xAndy, ThetaFromYTilde_NPlus2(intermediate));
    }

    return nPlus2Cache[xAndy];
}

const std::vector<MatrixTerm_Intermediate>& InteractionTermsFromXY(
        const std::string& xAndy) {
    if (intermediateCache.count(xAndy) == 0) {
        std::string x(xAndy.begin(), xAndy.begin() + xAndy.size()/2);
        std::string y(xAndy.begin() + xAndy.size()/2, xAndy.end());
        std::vector<char> uFromX(UFromX(x));
        std::vector<MatrixTerm_Intermediate> terms(YTildeFromY(y));
        for (auto& term : terms) {
            if (term.uPlus.size() < uFromX.size()/2) {
                term.uPlus.resize(uFromX.size()/2, 0);
                term.uMinus.resize(uFromX.size()/2, 0);
                term.yTilde.resize(uFromX.size()/2, 0);
            }
            for (std::size_t i = 0; i < term.uPlus.size(); ++i) {
                term.uPlus[i] += uFromX[i];
                term.uMinus[i] += uFromX[term.uPlus.size() + i];
                // term.coeff *= std::pow(std::sqrt(2), term.uPlus[i] + term.uMinus[i]);
            }
        }
        intermediateCache.emplace(xAndy, std::move(terms));
    }

    return intermediateCache[xAndy];
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
// isn't, so we have to subtract that off; that is, this computes Fbar, not F
std::string ExtractXY(const Mono& extractFromThis) {
    auto n = extractFromThis.NParticles();
    std::string xAndy(2*n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        xAndy[i] = extractFromThis.Pm(i) - 1;
        xAndy[n + i] = extractFromThis.Pt(i);
    }
    return xAndy;
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
        termsFromThisYTerm.back().coeff = yTerm.coeff;
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

        ret.back().coeff = coeff;
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

// the coefficient of a YTildeTerm, i.e. everything that's not a u or yTilde
coeff_class YTildeCoefficient(const char a, const char l, 
		const std::string& nAndm) {
    //coeff_class ret = ExactBinomial(a, l);
    coeff_class ret = Multinomial::Choose(2, a, {static_cast<char>(a-l),l});
    ret *= Multinomial::Lookup(nAndm.size()-1, nAndm);
    if ((a-l) % 2 == 1) ret = -ret;
    return ret;
}

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
        term.yTilde.resize(term.yTilde.size()-1); // so it can become "cosines"
        ret.emplace_back(term.coeff,
                std::move(term.uPlus), 
                std::move(term.uMinus),
                std::move(sines),
                std::move(term.yTilde) );
    }

    return ret;
}

// this is the same as above for all but the last two coordinates: for these,
// y(n-1)' -> sin(theta') and y(n)' -> cos(theta')
std::vector<MatrixTerm_Final> ThetaFromYTilde_NPlus2(
        std::vector<MatrixTerm_Intermediate>& intermediateTerms) {
    std::vector<MatrixTerm_Final> ret;

    for (auto& term : intermediateTerms) {
        // sine[i] appears in all yTilde[j] with j > i (strictly greater)
        std::vector<char> sines(term.yTilde.size()-2, 0);
        for (auto i = 0u; i < sines.size()-1; ++i) {
            for (auto j = i+1; j < term.yTilde.size(); ++j) {
                sines[i] += term.yTilde[j];
            }
        }

        // the two theta' functions (yTilde becomes the cosines)
        sines.back() = term.yTilde[term.yTilde.size()-2];
        term.yTilde[term.yTilde.size()-3] = term.yTilde[term.yTilde.size()-1];
        term.yTilde.resize(term.yTilde.size()-2);

        ret.emplace_back(term.coeff,
                std::move(term.uPlus), 
                std::move(term.uMinus),
                std::move(sines),
                std::move(term.yTilde) );
    }

    return ret;
}

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
            ret.emplace_back(term1.coeff * term2.coeff,
                    AddVectors(term1.uPlus, term2.uPlus),
                    AddVectors(term1.uMinus, term2.uMinus),
                    AddVectors(term1.sinTheta, term2.sinTheta),
                    AddVectors(term1.cosTheta, term2.cosTheta));
        }
    }
    //std::cout << "Combined Fs has " << ret.size() << " terms." << std::endl;
    return ret;
}

std::vector<InteractionTerm_Step2> CombineInteractionFs(
        const std::vector<MatrixTerm_Intermediate>& F1, 
        const std::vector<MatrixTerm_Intermediate>& F2) {
    std::vector<InteractionTerm_Step2> output;
    for (const auto& f1 : F1) {
        for (const auto& f2 : F2) {
            output.push_back(CombineInteractionFs_OneTerm(f1, f2));
        }
    }
    return output;
}

InteractionTerm_Step2 CombineInteractionFs_OneTerm(
        const MatrixTerm_Intermediate& f1, const MatrixTerm_Intermediate& f2) {
    InteractionTerm_Step2 output;
    output.coeff = f1.coeff * f2.coeff;

    output.u.resize(f1.uPlus.size() + f1.uMinus.size() + 2);
    for (std::size_t i = 0; i < f1.uPlus.size()-1; ++i) {
        output.u[2*i] = f1.uPlus[i] + f2.uPlus[i];
        output.u[2*i + 1] = f1.uMinus[i] + f2.uMinus[i];
    }
    output.u[output.u.size() - 4] = f1.uPlus.back();
    output.u[output.u.size() - 3] = f1.uMinus.back();
    output.u[output.u.size() - 2] = f2.uPlus.back();
    output.u[output.u.size() - 1] = f2.uMinus.back();

    output.theta.resize(f1.yTilde.size()-2 + f2.yTilde.size()-2, 0);
    // sine[i] appears in all yTilde[j] with j > i (strictly greater)
    for (auto i = 0u; i < f1.yTilde.size()-2; ++i) {
        for (auto j = i+1; j < f1.yTilde.size()-1; ++j) {
            output.theta[2*i] += f1.yTilde[j] + f2.yTilde[j];
        }
        output.theta[2*i + 1] = f1.yTilde[i] + f2.yTilde[i];
    }

    output.r[0] = 0;
    for (std::size_t i = 0; i < f1.yTilde.size()-1; ++i) {
        output.r[0] += f1.yTilde[i] + f2.yTilde[i];
    }
    output.r[1] = f1.yTilde.back();
    output.r[2] = f2.yTilde.back();
    return output;
}

// WARNING: if type == MAT_MASS this changes the MatrixTerm_Final vector by
// by changing each term's uMinus entries by -2
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
    //std::cout << "Returning FinalResult = " << totalFromIntegrals << std::endl;
    return totalFromIntegrals;
}

// the last two components of the uPlus and uMinus in combinedFs are u'+ & u'-;
// the last components of sinTheta and cosTheta in combinedFs are for theta'
coeff_class FinalResult_NPlus2(const std::vector<MatrixTerm_Final>& combinedFs){
    coeff_class totalFromIntegrals = 0;
    for (const auto& term : combinedFs) {
        totalFromIntegrals += DoAllIntegrals_NPlus2(term);
    }

    return totalFromIntegrals;
}

// do all of the integrals which are possible before mu discretization, and
// return a list of {value, {r exponents}} objects
//
// WARNING: this changes combinedFs so that it can't be reused
std::vector<InteractionTerm_Output> InteractionOutput(
        std::vector<InteractionTerm_Step2>& combinedFs, 
        const MATRIX_TYPE type, const coeff_class prefactor) {
    std::vector<InteractionTerm_Output> output;
    for (auto& combinedF : combinedFs) {
        output.emplace_back(prefactor*DoAllIntegrals(combinedF, type), 
                std::move(combinedF.r) );
    }
    return output;
}

// prefactors -----------------------------------------------------------------

coeff_class Prefactor(const Mono& A, const Mono&, const MATRIX_TYPE type) {
    if (type == MAT_INNER) {
        return InnerProductPrefactor(A.NParticles());
    } else if (type == MAT_MASS) {
        return MassMatrixPrefactor(A.NParticles());
    } else if (type == MAT_INTER_SAME_N) {
        return InteractionMatrixPrefactor(A.NParticles());
    } else if (type == MAT_INTER_N_PLUS_2) {
        return NPlus2MatrixPrefactor(A.NParticles());
    }
    std::cerr << "Error: prefactor type not recognized." << std::endl;
    return 0;
}

namespace {
    std::unordered_map<char, coeff_class> ipPrefactorCache;
    std::unordered_map<char, coeff_class> sameNPrefactorCache;
    std::unordered_map<char, coeff_class> nPlus2PrefactorCache;
} // anonymous namespace

// this follows (2.2) in Matrix Formulas.pdf
coeff_class InnerProductPrefactor(const char n) {
    if (ipPrefactorCache.count(n) == 0) {
        coeff_class denominator = std::tgamma(n+1); // tgamma = "true" gamma fcn
        denominator *= std::pow(16, n-1);
        denominator *= std::pow(M_PI, 2*n-3);
        //std::cout << "PREFACTOR: " << 1/denominator << std::endl;
        ipPrefactorCache.emplace(n, 1/denominator);
    }

    return ipPrefactorCache[n];
}

// this follows (2.3) in Matrix Formulas.pdf
coeff_class MassMatrixPrefactor(const char n) {
    // return n*InnerProductPrefactor(n); // if we're permuting M^2, remove n
    return InnerProductPrefactor(n);
}

coeff_class InteractionMatrixPrefactor(const char n) {
    if (sameNPrefactorCache.count(n) == 0) {
        coeff_class denominator = std::tgamma(n-1);
        denominator *= std::pow(M_PI*M_PI, n-1);
        denominator *= 4*std::pow(16, n);
        sameNPrefactorCache.emplace(n, 1/denominator);
    }

    return sameNPrefactorCache[n];
}

coeff_class NPlus2MatrixPrefactor(const char n) {
    if (nPlus2PrefactorCache.count(n) == 0) {
        coeff_class denominator = std::tgamma(n);
        denominator *= 6;
        denominator *= std::pow(M_PI, 2*n);
        denominator *= std::pow(16, n+1);
        nPlus2PrefactorCache.emplace(n, 1/denominator);
    }

    return nPlus2PrefactorCache[n];
}

// integrals ------------------------------------------------------------------

// do all the integrals for a direct matrix computation
coeff_class DoAllIntegrals(const MatrixTerm_Final& term) {
    std::size_t n = term.uPlus.size() + 1;
    coeff_class output = term.coeff;

    // do the u integrals first
    for (auto i = 0u; i < n-1; ++i) {
        output *= UIntegral(term.uPlus[i] + 3, 5*(n - (i+1)) - 2 + term.uMinus[i]);
    }

    // now the theta integrals; sineTheta.size() = cosTheta.size() = n-2.
    // All but the last one are short, while the last one is long
    //
    // these have constant terms which differ from Nikhil's because his i
    // starts at 1 instead of 0, so I use (i+1) instead
    if (n >= 3) {
        for (auto i = 0u; i < n-3; ++i) {
            output *= ThetaIntegral_Short(n - (i+1) - 2 + term.sinTheta[i],
                            term.cosTheta[i] );
        }
        output *= ThetaIntegral_Long(term.sinTheta[n-3], term.cosTheta[n-3]);
    }
    // std::cout << term.coeff << " * {" << term.uPlus << ", " << term.uMinus
        // << ", " << term.sinTheta << ", " << term.cosTheta << "} -> " << output 
        // << std::endl;
    return output;
}

// do all the integrals for an interaction matrix computation
//
// WARNING: this changes the exponents in term, rendering it non-reusable
//
// FIXME: may as well hash the results of this
coeff_class DoAllIntegrals(InteractionTerm_Step2& term, const MATRIX_TYPE type){
    std::cout << "DoAllIntegrals(" << term.u << ", " << term.theta << ")"
        << std::endl;
    auto n = term.u.size()/2;
    // using Nikhil's conventions, adjust exponents before doing integrals
    if (type == MAT_INTER_SAME_N) {
        for (std::size_t i = 0; i < n-2; ++i) {
            term.u[2*i] += 3;
            term.u[2*i + 1] += 5*(n-i) - 3;
        }
        term.u[term.u.size()-1] += 1;
        term.u[term.u.size()-2] += 1;
        term.u[term.u.size()-3] += 1;
        term.u[term.u.size()-4] += 1;

        for (std::size_t k = 0; k+4 < n; ++k) {
            term.theta[2*k] += n-k-3;
        }

        term.r[0] += n-3;
        term.r[1] -= 1;
        term.r[2] -= 1;
    } else if (type == MAT_INTER_N_PLUS_2) {
    }

    // actually do the integrals below
    coeff_class product = 1;
    for (std::size_t i = 0; i < n; ++i) {
        product *= UIntegral(term.u[2*i], term.u[2*i + 1]);
    }
    for (std::size_t k = 0; k < n-3; ++k) {
        product *= ThetaIntegral_Short(term.theta[2*k], term.theta[2*k + 1]);
    }
    return product;
}

// do all the integrals for an n+2 interaction computation
coeff_class DoAllIntegrals_NPlus2(const MatrixTerm_Final& term) {
    std::size_t n = term.uPlus.size() - 1;
    coeff_class output = term.coeff;

    // do the non-primed u integrals first
    for (auto i = 0u; i < n-1; ++i) {
        output *= UIntegral(term.uPlus[i] + 3, 5*(n - i) - 7 + term.uMinus[i]);
    }

    // next the two primed u integrals (FIXME: check additions!!)
    output *= UIntegral(term.uPlus[n-1] + 1, term.uMinus[n-1] + 1);
    output *= UIntegral(term.uPlus[n] + 1, term.uMinus[n] + 4);

    // now the theta integrals; there are n-2 "normal" theta integrals, followed
    // by one primed one. The primed and the last normal are long (I think?)
    //
    // these have constant terms which differ from Zuhair's because his i
    // starts at 1 instead of 0, so I use (i+1) instead
    if (n >= 3) {
        for (auto i = 0u; i < n-3; ++i) {
            output *= ThetaIntegral_Short(n - (i+1) - 2 + term.sinTheta[i],
                            term.cosTheta[i] );
        }
        output *= ThetaIntegral_Long(term.sinTheta[n-3], term.cosTheta[n-3]);
    }
    // I think this one (the primed one) is guaranteed to be there
    output *= ThetaIntegral_Long(term.sinTheta[n-2], term.cosTheta[n-2]);

    // std::cout << term.coeff << " * {" << term.uPlus << ", " << term.uMinus
        // << ", " << term.sinTheta << ", " << term.cosTheta << "} -> " << output 
        // << std::endl;
    return output;
}

// this is the integral over the "u" variables (u \el [-1,1]; it uses a 2F1
// identity to turn 2F1(a, b, c, -1) -> 2^(-a)*2F1(a, c-b, c, 1/2)
//
// follows the conventions of Zuhair's 5.34; a is the exponent of u_i+ and
// b is the exponent of u_i-
//
// this can also be thought of as the integral over z, where a is the exponent
// of sqrt(z) and b is the exponent of sqrt(1 - z); in this case z \el [0,1].
builtin_class UIntegral(const builtin_class a, const builtin_class b) {
    //std::cout << "UIntegral(" << a << ", " << b << ")" << std::endl;
    std::array<builtin_class,2> abArray{{a,b}};
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
builtin_class ThetaIntegral_Short(const builtin_class a, const builtin_class b) {
    if (static_cast<int>(b) % 2 == 1) return 0;
    std::array<builtin_class,2> abArray{{a,b}};
    if (b < a) std::swap(abArray[0], abArray[1]);
    if (thetaCache.count(abArray) == 1) return thetaCache.at(abArray);

    builtin_class ret = std::exp(std::lgamma((1+a)/2) + std::lgamma((1+b)/2) 
                    - std::lgamma((2 + a + b)/2) );
    thetaCache.emplace(abArray, ret);
    return ret;
}

// this is the integral over the "theta" veriables from 0 to 2pi; it implements 
// Zuhair's 5.36, where a is the exponent of sin(theta) and b is the exponent of
// cos(theta).
builtin_class ThetaIntegral_Long(const builtin_class a, const builtin_class b) {
    if (static_cast<int>(a+b) % 2 == 1) return 0;
    return 2*ThetaIntegral_Short(a, b);
}

// Integral of r -- the answer depends on alpha, so I think this can not be
// completed until the discretization step
//
// Actually, the alpha dependence can be recorded: it depends on a via 
// arguments to the hypergeometric functions, and if 0 < alpha < 1 then we get
// the alpha > 1 answer with alpha -> 1/alpha, except that a factor of 
// alpha^(-1 - a) disappears
//
// The answer given assumes alpha > 1, so it must start there or be inverted.
// The other answer is probably actually cleaner, since it has alpha^2 in the
// hypergeometrics instead of 1/alpha^2.
/*builtin_class RIntegral(const builtin_class a, builtin_class alpha) {
    if (alpha < 1) {
        rCache.emplace(alpha, RIntegral(a, 1/alpha)*std::pow(alpha, -a-1));
    }
    if (rCache.count(a) == 0) {
        value = std::sqrt(M_PI) * std::tgamma((a+1)/2) / 4;
        value /= (alpha*alpha - 1);
        // if (alpha > 1) value *= std::pow(alpha, -a - 1); // OBVIOUSLY WON'T WORK

        builtin_class hypergeos = 2 * alpha*alpha *
            gsl_sf_hyperg_2F1(-0.5, (a+1)/2, a/2+1, 1/(alpha*alpha))
            / std::tgamma(a/2 + 1);
        hypergeos -= gsl_sf_hyperg_2F1(0.5, (a+1)/2, a/2+2, 1/(alpha*alpha))
            / std::tgamma(a/2 + 2);
        value *= hypergeos;
        rCache.emplace(aalpha, value);
    }

    return rCache.at(a);
}*/

} // namespace MatrixInternal
