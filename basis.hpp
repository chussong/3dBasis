#ifndef BASIS_HPP
#define BASIS_HPP

#include <vector>
#include <sstream>

#include "constants.hpp"
#include "construction.hpp"
#include "mono.hpp"
#include "poly.hpp"

// a basis is anything which claims that it can express monomials and
// polynomials as vectors on itself. It can contain monos or polys as needed.
// * Can be accessed with [] like an array and iterated through with begin()
// and end().
// * Has a stream operator, which outputs all the basis vectors on new lines.
// * Most important function is ExpressPoly, which takes a polynomial and
// expresses it as a vector in this basis's space.
template<class T> class Basis;
template<class T> std::ostream& operator<<(std::ostream& os, const Basis<T>& out);

template<class T>
class Basis {
    std::vector<T> basisVectors;

    public:
        explicit Basis(const std::vector<T>& basisVectors): basisVectors(basisVectors) {}
        Basis(const int numP, int degree, const Arguments& args);
        //explicit Basis(const Basis&) = default;
        //Basis(Basis&&) = default;

        unsigned int FindInBasis(const std::vector<int>& pm,
                        const std::vector<int>& pt) const;
        unsigned int FindInBasis(const T& wildVector) const;

        void DeleteOdd();
        void DeleteEven();

        const T& operator[](size_t i) const { return basisVectors[i]; }
        friend std::ostream& operator<<<T>(std::ostream& os, const Basis<T>& out);

        std::size_t size() const { return basisVectors.size(); }
        typename std::vector<T>::const_iterator		begin() const noexcept
                        { return basisVectors.begin(); }
        typename std::vector<T>::iterator			begin() noexcept
                        { return basisVectors.begin(); }
        typename std::vector<T>::const_iterator		end()	const noexcept
                        { return basisVectors.end(); }
        typename std::vector<T>::iterator			end()	noexcept
                        { return basisVectors.end(); }

        // Triplet ExpressMono(const Mono& toExpress, const int column,
                        // const int rowOffset) const;
        // std::list<Triplet> ExpressPoly(const Poly& toExpress, 
                        // const int column, const int rowOffset) const;
        DVector DenseExpressMono(const Mono& toExpress) const;
        DVector DenseExpressPoly(const Poly& toExpress) const;

};

// Contains two bases and intelligently decides which one to use for various
// requests. Because this is in two pieces, it can not be iterated through.
template<class T>
class splitBasis {
    Basis<T> evenBasis;
    Basis<T> oddBasis;

    public:
        splitBasis(const int numP, const int degree, const Arguments& args);

        std::pair<unsigned int, Basis<T>*> FindInBasis(const std::vector<int>& pm,
                        const std::vector<int>& pt);
        std::pair<unsigned int, Basis<T>*> FindInBasis(const Mono& wildMono);

        Basis<T>& OddBasis() { return oddBasis; }
        const Basis<T>& OddBasis() const { return oddBasis; }
        Basis<T>& EvenBasis() { return evenBasis; }
        const Basis<T>& EvenBasis() const { return evenBasis; }

        // std::list<Triplet> ExpressPoly(const Poly& toExpress, const int column,
                        // const int row) const;

        static bool IsOdd (const Mono& toTest);
        static bool IsEven(const Mono& toTest);
        static std::pair<Poly, Poly> OddEvenSplit  (const Poly& toSplit);
};

// state generation -----------------------------------------------------------

bool EoMAllowed(const std::vector<particle>& cfg);

template<class T>
inline Basis<T>::Basis(const int, int, const Arguments&) {
    std::cerr << "Error: ordered to construct a basis by degree for an "
        << "underlying type which has not been specialized. Please construct "
        << "with a different type or write a specialization for this type."
        << std::endl;
}

template<>
inline Basis<Mono>::Basis(const int numP, int degree, const Arguments& args) {
    // 1: generate all possibilities for P_-
    // 2: identify nodes in P_-, generate possible distributions of P_\perp to
    // 		the nodes and then P_\perp within each node, adding each at top level
    // 3: identify nodes in (P_- and P_\perp), repeating step 2 for P_+
    // 4: take list from step 3 and create a mono from each entry, then store
    // 		the list of these as basisMonos
    // NOTE: this does not make any use of the symmetries between the components
    // 		in constructing the basis; for instance, the states where every
    // 		pm = 0 are followed by a copy of the earlier parts of the basis

    int options = args.options;
    OStream& console = *args.console;
    const bool debug = ((options & OPT_DEBUG) != 0);
    if(debug) console << "***Generating basis at N=" << numP << ", D="
            << degree << "***" << endl;

    // subtract off the required Dirichlet derivatives
    degree -= numP;

    if(degree < 0){
        console << "Error: there are no Dirichlet states with degree "
                << "<= the number of particles." << endl;
        return;
    }

    // start with all possible configurations of P_-
    if(degree == 0){
        std::vector<particle> onlyMono;
        onlyMono.resize(numP);
        for(auto& part : onlyMono){
            part.pm = 1;
        }
        basisVectors.emplace_back(onlyMono);
        return;
    }
    std::vector<std::vector<int>> minus;
    if(options & OPT_ALLMINUS){
        minus = GetStatesAtDegree<int>(numP, degree);
    } else {
        minus = GetStatesUpToDegree<int>(numP, degree);
    }
    std::vector<std::vector<particle>> particleCfgs;
    for(auto& minusCfg : minus){
        if(debug) console << "Here's a configuration of minuses: "
            << minusCfg << endl;
        std::vector<particle> newCfg(minusCfg.size());
        for(auto i = 0u; i < newCfg.size(); ++i) newCfg[i].pm = minusCfg[i];
        particleCfgs.push_back(newCfg);
    }

    // for each P_- configuration, generate all possible P_\perp configs
    std::vector<int> nodes;
    std::vector<std::vector<particle>> newCfgs;
    for(auto& configuration : particleCfgs){
        if(options & OPT_ALLMINUS){
            newCfgs.push_back(configuration);
            continue;
        }
        nodes = IdentifyNodes(configuration);
        int remainingEnergy = degree;
        for(auto& part : configuration) remainingEnergy -= part.pm;
        std::vector<std::vector<int>> perp(CfgsFromNodes(remainingEnergy, nodes,
                                                                        true));
        if(debug) console << "COMBINING SUBCONFIGS OF " << configuration 
            << endl;
        for(auto& newCfg : CombinedCfgs(configuration, perp, 2)){
            if(debug) console << "NEW CONFIGURATION: " << newCfg << endl;
            newCfgs.push_back(newCfg);
        }
    }

    // finally, add in the required P_- from being Dirichlet
    particleCfgs.clear();
    for(auto& cfg : newCfgs){
        for(auto& part : cfg) ++part.pm;
        particleCfgs.push_back(cfg);
    }

    //console << "Tick." << endl;
    for(auto& cfg : particleCfgs){
        basisVectors.emplace_back(cfg);
        //console << cfg << endl;
    }
}

template<class T>
inline unsigned int Basis<T>::FindInBasis(const std::vector<int>&,
		const std::vector<int>&) const{
    std::cerr << "Unspecialized Basis::FindInBasis(vectors) should not be called. "
        << "How was this basis object created in the first place?" << std::endl;
    return -1u;
}

/*
template<class T>
inline Triplet Basis<T>::ExpressMono(const Mono&, const int, const int) const{
    std::cout << "Unspecialized Basis::ExpressMono should never be called. "
        << "How was a basis object created in the first place?" << endl;
    return Triplet(-1, -1, coeff_class(0));
}

template<class T>
inline std::list<Triplet> Basis<T>::ExpressPoly(const Poly&, const int, const int) const{
    std::cout << "basis::ExpressMono should never be called. "
        << "How was a basis object created in the first place?" << endl;
    return {{Triplet(-1, -1, coeff_class(0))}};
}
*/

template<class T>
inline DVector Basis<T>::DenseExpressMono(const Mono&) const {
    std::cerr << "Basis<T>::DenseExpressMono should never be called. "
        << "How was a Basis<T> object created in the first place?" << std::endl;
    return DVector(0);
}

template<class T>
inline DVector Basis<T>::DenseExpressPoly(const Poly& toExpress) const {
    std::cerr << "Basis<T>::DenseExpressPoly should never be called. "
        << "How was a Basis<T> object created in the first place?" << std::endl;
    return DVector(toExpress.size());
}

template<class T>
inline unsigned int Basis<T>::FindInBasis(const T& wildVector) const{
    for(auto i = 0u; i < basisVectors.size(); ++i){
        if(basisVectors[i] == wildVector) return i;
    }
    std::cerr << "Warning! Failed to find the following vector in our basis: "
        << wildVector << std::endl;
    return -1u;
}

template<>
inline unsigned int Basis<Mono>::FindInBasis(const std::vector<int>& pm, 
		const std::vector<int>& pt) const{
    return FindInBasis(Mono(pm, pt));
}

template<>
inline void Basis<Mono>::DeleteOdd(){
    basisVectors.erase(std::remove_if(basisVectors.begin(), basisVectors.end(), 
                            splitBasis<Mono>::IsOdd), basisVectors.end());
}

template<>
inline void Basis<Mono>::DeleteEven(){
    basisVectors.erase(std::remove_if(basisVectors.begin(), basisVectors.end(),
                            splitBasis<Mono>::IsEven), basisVectors.end());
}

template<class T>
inline std::ostream& operator<<(std::ostream& os, const Basis<T>& out){
    // std::cout << "Basis<Mono> goin' out" << std::endl;
    if (out.size() == 0) return os << "{ }";
    os << "{ ";
    for (std::size_t i = 0; i < out.size()-1; ++i) {
        os << out[i].HumanReadable() << ", ";
    }
    return os << out[out.size()-1].HumanReadable() << " }";
}

template<class T>
inline std::string MathematicaOutput(const Basis<T>& out) {
    if(out.size() == 0) return "{ }";
    std::stringstream ss;
    ss << "{ ";
    for (std::size_t i = 0; i < out.size()-1; ++i) {
        ss << MathematicaOutput(out[i]) << ", ";
    }
    ss << MathematicaOutput(out[out.size()-1]) << " }";
    return ss.str();
}

/*
template<>
inline Triplet Basis<Mono>::ExpressMono(const Mono& toExpress, const int column,
		const int rowOffset) const{
    for(auto i = 0u; i < basisVectors.size(); ++i){
        if(toExpress == basisVectors[i])
            return Triplet(rowOffset+i, column, 
                            toExpress.Coeff()/basisVectors[i].Coeff());
    }
    std::cerr << "Error: tried to express the monomial " << toExpress
    << " on the given basis but was not able to identify it." << std::endl;
    return Triplet(-1, -1, toExpress.Coeff());
}

template<>
inline std::list<Triplet> Basis<Mono>::ExpressPoly(const Poly& toExpress, 
		const int column, const int rowOffset) const{
    std::list<Triplet> ret;
    if(toExpress.size() == 0) return ret;
    unsigned int hits = 0u;
    for(auto i = 0u; i < basisVectors.size(); ++i){
        for(auto& term : toExpress){
            // if(std::abs(term.Coeff()) < EPSILON){
                // ++zeros;
                // std::cout << "Zero get: " << term.Coeff() << ". Now have "
                        // << zeros << "." << std::endl;
                // continue;
            // } it should not actually be possible for a Poly to have coeff = 0
            if(term == basisVectors[i]){
                ret.emplace_front(rowOffset+i, column, 
                                term.Coeff()/basisVectors[i].Coeff());
                ++hits;
                break;
            }
        }
        if(hits == toExpress.size()) break; // right??
    }
    if(hits < toExpress.size()){
        std::cerr << "Error: tried to express the polynomial " << toExpress
        << " on the given basis but was only able to identify " << hits
        << " of " << toExpress.size() << " terms. Here's the basis:" << std::endl;
        std::cerr << *this << std::endl;
    }
    //std::cout << "Expressed " << toExpress << " as the following triplets:" << std::endl;
    //for(auto& trip : ret) std::cout << trip << std::endl;
    return ret;
}
*/

template<>
inline DVector Basis<Mono>::DenseExpressMono(const Mono& toExpress) const {
    DVector output = DVector::Zero(this->size());
    for (std::size_t i = 0; i < this->size(); ++i) {
        if ((*this)[i] == toExpress) {
            output(i) = toExpress.Coeff();
            return output;
        }
    }
    std::cerr << "Error: attempted to express " << toExpress << " on the basis "
            << *this << " but was not able to." << std::endl;
    return output;
}

template<>
inline DVector Basis<Mono>::DenseExpressPoly(const Poly& toExpress) const {
    DVector output = DVector::Zero(this->size());
    for (const auto& term : toExpress) output += DenseExpressMono(term);
    return output;
}

/*
// This does not attempt to find a linear combination of known polynomials
// which would reproduce toExpress. Obviously it would be more correct if it
// did attempt to do so.
template<>
inline Triplet Basis<Poly>::ExpressMono(const Mono& toExpress, const int column,
		const int rowOffset) const{
    Poly polyForm(toExpress);
    for(auto i = 0u; i < basisVectors.size(); ++i){
        if(polyForm == basisVectors[i])
                return Triplet(rowOffset+i, column, 
                                toExpress.Coeff()/basisVectors[i][0].Coeff());
    }
    std::cerr << "Error: tried to express the monomial " << toExpress
    << " on the given basis but was not able to identify it." << std::endl;
    return Triplet(-1, -1, toExpress.Coeff());
}

// Again, no linear combinations are checked -- we just look to see if we have
// the poly as written. This is obviously not ideal. It also means we can only
// get exactly 1*the polynomial, which actually probably breaks stuff.
template<>
inline std::list<Triplet> Basis<Poly>::ExpressPoly(const Poly& toExpress, 
		const int column, const int rowOffset) const{
    for(auto i = 0u; i < basisVectors.size(); ++i){
            if(toExpress == basisVectors[i])
                    return {{Triplet(rowOffset+i, column, 1)}};
    }
    std::cerr << "Error: tried to express the plynomial " << toExpress
    << " on the given basis but was not able to identify it." << std::endl;
    return {{Triplet(-1, -1, 0)}};
}
*/

template<class T>
inline bool splitBasis<T>::IsOdd(const Mono& toTest){
    return toTest.TotalPt()%2 == 1;
}

template<class T>
inline bool splitBasis<T>::IsEven(const Mono& toTest){
    return !IsOdd(toTest);
}

// this could definitely be done more intelligently if speed were important;
// even at order 0, we could just generate the basis once and copy it rather
// than making exactly the same thing twice.
template<class T>
inline splitBasis<T>::splitBasis(const int numP, const int degree, 
        const Arguments& args): evenBasis(numP, degree, args), 
                                oddBasis(evenBasis) {
    //std::cout << evenBasis << "----->" << std::endl;
    evenBasis.DeleteOdd();
    //std::cout << evenBasis << std::endl;
    oddBasis.DeleteEven();
    //std::cout << oddBasis << std::endl;
}

/*
template<class T>
inline std::list<Triplet> splitBasis<T>::ExpressPoly(const Poly& toExpress, 
		const int column, const int rowOffset) const{
    std::list<Triplet> ret;
    for(auto& term : toExpress){
        if(IsOdd(term)){
            ret.push_front(oddBasis.ExpressMono(term, column, rowOffset));
        } else {
            ret.push_front(evenBasis.ExpressMono(term, column, rowOffset));
        }
    }
    return ret;
}
*/

// Priority for sorting monomials within a basis. The point of this is to get
// the orthogonalization that we want from Gram-Schmidt, so it favors states 
// with lower degree, and then states with lower total p_\perp as a tiebreaker;
// within these catergories, states with larger p_- and p_\perp hierarchies are
// favored
inline bool SortPriority(const Mono& A, const Mono& B) {
    if (A.Degree() != B.Degree()) {
        return A.Degree() < B.Degree();
    } else if (A.TotalPt() != B.TotalPt()) {
        return A.TotalPt() < B.TotalPt();
    } else if (A.MaxPm() != B.MaxPm()) {
        // within a given degree, favor states that are tall rather than wide
        return A.MaxPm() > B.MaxPm();
    } else {
        return A.MaxPt() > B.MaxPt();
    }
}

//bool SortPriority(const Poly& A, const Poly& B);

// combine multiple bases into a single one
template<class T>
Basis<T> CombineBases(const std::vector<Basis<T>>& oldBases) {
    std::vector<T> newBasisVectors;
    for (auto& oldBasis : oldBases) {
        for (auto& oldBasisVector : oldBasis) {
            newBasisVectors.push_back(oldBasisVector);
        }
    }
    std::sort(newBasisVectors.begin(), newBasisVectors.end(), SortPriority);
    return Basis<T>(newBasisVectors);
}

// return the minimal basis of monomials needed to express the given polynomials
inline Basis<Mono> MinimalBasis(const std::vector<Poly>& polynomials) {
    Poly combinedPoly;
    for (const auto& poly : polynomials) {
        combinedPoly += poly;
    }
    std::vector<Mono> allUsedMonos(combinedPoly.size());
    for (auto i = 0u; i < combinedPoly.size(); ++i) {
        allUsedMonos[i] = combinedPoly[i];
        allUsedMonos[i].Coeff() = 1;
    }
    return Basis<Mono>(allUsedMonos);
}

// get an element from a vector of multiple bases treated like a single basis
template<class T>
const T& Get(const std::vector<Basis<T>>& multipleBases, size_t index){
    for(const auto& basis : multipleBases){
        if(index < basis.size()){
            return basis[index];
        } else {
            index -= basis.size();
        }
    }
    throw std::runtime_error("out of range error in Get(vector<Basis>)");
}

// Take a basis and look at every element; if the element has zero norm (within
// floating point tolerance), delete it, otherwise normalize it.
template<class T>
void Normalize(Basis<T>& toNormalize) {
    // if vector's norm is zero, skip it; otherwise, normalize
    std::vector<T> newBasisVectors;
    for (const T& basisVector : toNormalize) {
        if (basisVector.IsNull()) {
            //std::cout << basisVector.HumanReadable() << " judged to be a null "
                    //<< "state." << std::endl;
            continue;
        }
        builtin_class norm = InnerFock(basisVector, basisVector);
        norm = std::sqrt(norm);
        // coeff_class norm = T::InnerProduct(basisVector, basisVector, cache, 
                        // kCache);
        /*std::cout << "Norm of " << basisVector.HumanReadable() << ": " 
                << norm << std::endl;
        std::cout << "Meanwhile, its Kt is: " << basisVector.K2(0.5) << "."
                << std::endl;*/
        // using norm instead of std::abs(norm) because negative norms are zeros
        //basisVector /= std::sqrt(norm);
        newBasisVectors.push_back(basisVector);
        newBasisVectors.back() /= norm;
    }
    Basis<T> newBasis{newBasisVectors};
    std::swap(toNormalize, newBasis);

    // delete everything with a zero norm
    /*toNormalize.erase(std::remove_if(toNormalize.begin(), toNormalize.end(),
                            [](const T& vec){return vec.Coeff() == 0;}), 
                    toNormalize.end());*/
}

// take vector in monomial or polynomial basis, return it as polynomial
template<class T>
Poly VectorToPoly(const DVector& kernelVector, const Basis<T>& startBasis){
    Poly ret;
    if(static_cast<size_t>(kernelVector.rows()) != startBasis.size()){
        std::cerr << "Error: the given Q column has " << kernelVector.rows()
            << " rows, " << "but the given basis has " << startBasis.size() 
            << " monomials. These must be the same." << std::endl;
        return ret;
    }
    for(auto row = 0; row < kernelVector.rows(); ++row){
        if(std::abs<builtin_class>(kernelVector.coeff(row)) < EPSILON) continue;
        ret += kernelVector.coeff(row)*startBasis[row];
    }

    // if(ret.size() == 0) return ret;
    // coeff_class smallestCoeff = std::abs(ret[0].Coeff());
    // for(auto& term : ret) smallestCoeff = 
            // std::min(std::abs(term.Coeff()), smallestCoeff);
    // for(auto& term : ret) term /= smallestCoeff;
    return ret;
}

#endif
