#include "3dBasis.hpp"

int main(int argc, char* argv[]) {
	arguments args = ParseArguments(argc, argv);
	if (args.options & OPT_VERSION){
		std::cout << "This is 3dBasis version " << VERSION << ", released "
			<< RELEASE_DATE << ". The latest updates can always be found at "
			<< "https://github.com/chussong/3dBasis." << std::endl;
		return EXIT_SUCCESS;
	}

	if (args.degree == 0 || args.numP == 0){
		std::cerr << "Error: you must enter a number of particles, a degree, "
			<< "and a value for delta." << std::endl;
		return EXIT_FAILURE;
	}

	if (args.options & OPT_MULTINOMTEST) {
		Multinomial::Initialize(args.numP, args.degree);
		for (char n = 0; n <= args.degree; ++n) {
			std::cout << "n = " << std::to_string(n) << ": ";
			for (auto& mVector : Multinomial::GetMVectors(n)) {
				//std::cout << std::endl << Multinomial::MVectorOut(mVector) << ": ";
				std::cout << Multinomial::Lookup(mVector) << ", ";
			}
			std::cout << std::endl;
		}
		return EXIT_SUCCESS;
	}

	if (std::abs(args.delta) < EPSILON) args.delta = 0.5;

	//if(args.options & OPT_IPTEST){
		return InnerProductTest(args);
	//}
}

int InnerProductTest(const arguments& args){
	int numP = args.numP;
	int degree = args.degree + args.numP; // add required Dirichlet derivatives
	//coeff_class delta = args.delta;
	int options = args.options;

	//options = options | OPT_DEBUG;

	std::cout << "Beginning inner product test with N=" << numP << ", L="
		<< degree << "." << std::endl;
	
	//std::cout << "Testing gamma cache construction." << std::endl;
	GammaCache cache(numP, 2*degree, 2*(degree-numP));
	KVectorCache kCache(numP, 2*(degree-numP)); // maxPt might be half this?
	//std::cout << "A coefficient from it: " << cache.Middle(degree-2, 2, 1)
		//<< std::endl;

	std::vector<Basis<Mono>> allEvenBases;
	std::vector<Basis<Mono>> allOddBases;
	for(int deg = numP; deg <= degree; ++deg){
		splitBasis<Mono> degBasis(numP, deg, options);
		allEvenBases.push_back(degBasis.EvenBasis());
		allOddBases.push_back(degBasis.OddBasis());
		//std::cout << allEvenBases.back() << std::endl;
		//std::cout << allOddBases.back() << std::endl;
	}

	/*for(auto& basis : allEvenBases){
		for(auto& basisMono : basis){
			basisMono /= std::sqrt(Mono::InnerProduct(basisMono, basisMono,
						cache, kCache) );
		}
	}
	for(auto& basis : allOddBases){
		for(auto& basisMono : basis){
			basisMono /= std::sqrt(Mono::InnerProduct(basisMono, basisMono,
						cache, kCache) );
		}
	}*/

	std::cout << "EVEN STATE ORTHOGONALIZATION" << std::endl;
	Orthogonalize(allEvenBases, cache, kCache);
	std::cout << "REORDERED" << std::endl;
	std::swap(allEvenBases.front(), allEvenBases.back());
	Orthogonalize(allEvenBases, cache, kCache);

	std::cout << "ODD STATE ORTHOGONALIZATION" << std::endl;
	Orthogonalize(allOddBases, cache, kCache);
	std::cout << "REORDERED" << std::endl;
	std::swap(allOddBases.front(), allOddBases.back());
	Orthogonalize(allOddBases, cache, kCache);

	return EXIT_SUCCESS;
}

int Orthogonalize(const std::vector<Basis<Mono>>& inputBases,
					const GammaCache& cache, const KVectorCache& kCache){
	Timer timer;
	Basis<Mono> unifiedBasis = CombineBases(inputBases);
	Normalize(&unifiedBasis, cache, kCache);
	DMatrix gram = GramMatrix(unifiedBasis, cache, kCache);
	if(gram.rows() == 0) return 0;
	
	std::cout << "Gram matrix constructed in " << timer.TimeElapsedInWords()
		<< "." << std::endl;
	if(TotalSize(inputBases) <= 7){
		std::cout << gram << std::endl;
		//std::cout << "Compare to the non-unified one: " << std::endl 
			//<< GramMatrix(inputBases, cache, kCache) << std::endl;
	}

	/*timer.Start();
	DQRSolver solver;
	// it's possible there's a smarter way to set the threshold than this,
	// which could possibly take a long time? We multiply our custom epsilon by
	// the largest coefficient of gram's component-wise absolute value
	solver.setThreshold(EPSILON*gram.cwiseAbs().maxCoeff());
	//std::cout << "Matrix threshold set to " << EPSILON*gram.cwiseAbs().maxCoeff()
		//<< "." << std::endl;
	solver.compute(gram);
	DMatrix QMatrix = ExtractQMatrix(solver, gram.rows());
	std::cout << "Q matrix found in " << timer.TimeElapsedInWords() << "."
		<< std::endl;
	if(TotalSize(inputBases) <= 7){
		std::cout << QMatrix << std::endl;
	}
	std::cout << "Rank, i.e. number of independent operators: " << solver.rank() 
		<< std::endl;*/

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

	timer.Start();
	Basis<Mono> minimalBasis(MinimalBasis(orthogonalized));
	DMatrix massMatrix(MassMatrix(minimalBasis));

	std::cout << "Computed this mass matrix from the basis in " 
		<< timer.TimeElapsedInWords() << ", getting this matrix:" << std::endl;
	std::cout << massMatrix << std::endl;

	/*std::cout << "Here's the new gram matrix to confirm orthonormality:"
		<< std::endl;
	DMatrix gram2 = GramMatrix(Basis<Poly>(orthogonalized), cache, kCache);
	std::cout << gram2 << std::endl;*/

	//timer.Start();
	//std::vector<Poly> matrixGS = MatrixGramSchmidt(gram, unifiedBasis);
	//std::cout << "Matrix Gram-Schmidt performed in " 
		//<< timer.TimeElapsedInWords() << ", giving this:" << std::endl;
	//for(auto& p : matrixGS) std::cout << p << std::endl;

	return orthogonalized.size();
}

std::vector<Poly> GramSchmidt_WithMatrix(const std::vector<Basis<Mono>> input,
		const DMatrix& gramMatrix) {
	return GramSchmidt_WithMatrix(CombineBases(input), gramMatrix);
}

// this may be assuming that the incoming basis is already normalized (though
// of course not orthogonal)
std::vector<Poly> GramSchmidt_WithMatrix(const Basis<Mono> inputBasis, 
		const DMatrix& gramMatrix) {
	std::vector<DVector> vectorForms;
	//std::vector<coeff_class> inverseNorms;
	for (auto i = 0u; i < inputBasis.size(); ++i) {
		//DVector nextVector = DVector::Zero();
		//nextVector(i) = 1;
		DVector nextVector = DVector::Unit(inputBasis.size(), i);
		for (auto j = 0u; j < vectorForms.size(); ++j) {
			//nextVector -= gramMatrix(i, j)*vectorForms[j]*inverseNorms[j];
			nextVector -= GSProjection(nextVector, vectorForms[j], gramMatrix);
					//vectorForms[j]*inverseNorms[j], gramMatrix);
		}

		coeff_class norm = GSNorm(nextVector, gramMatrix);
		if (std::abs(norm) < EPSILON) continue;
		//nextVector /= std::sqrt(norm);
		//std::cout << "Projections of final vector onto old ones: " << std::endl;
		//for (auto& oldVec : vectorForms) {
			//std::cout << GSNorm(GSProjection(nextVector, oldVec, gramMatrix), gramMatrix)
				//<< std::endl;
		//}
		vectorForms.push_back(nextVector/std::sqrt(norm));
		//inverseNorms.push_back(1/norm);
	}

	std::vector<Poly> ret;
	for (auto& vec : vectorForms) ret.push_back(VectorToPoly(vec, inputBasis));
	return ret;
}

// projectOnto must be normalized already if you want the right answer from this
DVector GSProjection(const DVector& toProject, const DVector& projectOnto,
		const DMatrix& gramMatrix) {
	coeff_class scale = 0;
	for (auto i = 0u; i < toProject.size(); ++i) {
		if (toProject(i) == 0) continue;
		for (auto j = 0u; j < toProject.size(); ++j) {
			scale += toProject(i)*projectOnto(j)*gramMatrix(i,j);
		}
	}
	//std::cout << "Projection of \n" << toProject << " onto \n" << projectOnto
		//<< " is \n" << result << std::endl;
	return scale*projectOnto;
}

coeff_class GSNorm(const DVector& vector, const DMatrix& gramMatrix) {
	coeff_class norm = 0;
	for(auto i = 0u; i < vector.size(); ++i) {
		norm += vector(i)*vector(i)*gramMatrix(i,i);
		for (auto j = i+1; j < vector.size(); ++j) {
			norm += 2*vector(i)*vector(j)*gramMatrix(i,j);
		}
	}
	//std::cout << "This vector's norm is " << norm << ":" << std::endl << vector 
		//<< std::endl;
	return norm;
}

std::vector<Poly> GramSchmidt_MatrixOnly(const DMatrix& input, 
		const std::vector<Basis<Mono>>& inputBases) {
	std::vector<DVector> knownVectors;
	for(Eigen::Index col = 0; col < input.cols(); ++col){
		DVector newVector = input.col(col);
		for(const DVector& knownVector : knownVectors){
			newVector -= newVector.dot(knownVector)*knownVector;
		}
		if(newVector.norm() <= EPSILON) continue;
		newVector /= newVector.norm();
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

arguments ParseArguments(int argc, char* argv[]){
	std::vector<std::string> options;
	std::string arg;
	arguments ret;
	int j = 0;
	for(int i = 1; i < argc; ++i){
		arg = argv[i];
		if(arg.size() > 0){
			if(arg[0] == '-'){
				options.push_back(arg);
			} else {
				switch(j){
					case 0:
						ret.numP = ReadArg<int>(arg);
						break;
					case 1:
						ret.degree = ReadArg<int>(arg);
						break;
					case 2:
						ret.delta = ReadArg<coeff_class>(arg);
						break;
					default:
						std::cerr << "Error: at most two non-option arguments may "
							<< "be given." << std::endl;
						return ret;
				}
				++j;
			}
		}
	}
	if(j < 2) ret.numP = 0; // invalidate the input since it was insufficient
	ret.options = ParseOptions(options);
	return ret;
}

// -b solves using the non-split method
int ParseOptions(std::vector<std::string> options){
	int ret = 0;
	for(auto& opt : options){
		if(opt.compare(0, 2, "-v") == 0){
			ret = ret | OPT_VERSION;
			continue;
		}
		if(opt.compare(0, 2, "-d") == 0){
			ret = ret | OPT_DEBUG;
			ret = ret | OPT_OUTPUT;
			continue;
		}
		if(opt.compare(0, 2, "-o") == 0){
			ret = ret | OPT_OUTPUT;
			continue;
		}
		if(opt.compare(0, 2, "-i") == 0){
			ret = ret | OPT_IPTEST;
			continue;
		}
		if(opt.compare(0, 2, "-m") == 0){
			ret = ret | OPT_MULTINOMTEST;
			continue;
		}
		if(opt.compare(0, 2, "-M") == 0){
			ret = ret | OPT_ALLMINUS;
			continue;
		}
		if(opt.compare(0, 1, "-") == 0){
			std::cerr << "Warning: unrecognized option " << opt << " will be "
				<< "ignored." << std::endl;
			continue;
		}
	}
	return ret;
}

bool particle::operator==(const particle& other) const{
	return (pm == other.pm) && (pt == other.pt);
}

// note: triplets displayed (row, column, value) despite matrix's storage type
std::ostream& operator<<(std::ostream& os, const Triplet& out){
	return os << "(" << out.row() << "," << out.col() << "," << out.value()
		<< ")";
}

// this version is a 'loose' EoM compliance that only removes Pp which are on
// the same particle as a Pm
bool EoMAllowed(){
	std::cout << "Warning: EoMAllowed is deprecated." << std::endl;
	return true;
}

std::list<Triplet> ConvertToRows(const std::vector<Poly>& polyForms, 
		const Basis<Mono>& targetBasis, const Eigen::Index rowOffset){
	if(polyForms.size() == 0) return std::list<Triplet>();
	std::list<Triplet> ret = targetBasis.ExpressPoly(polyForms[0], 0,
			rowOffset);
	for(auto i = 1u; i < polyForms.size(); ++i){
		ret.splice(ret.end(), targetBasis.ExpressPoly(polyForms[i], i,
				rowOffset));
	}
	return ret;
}

Poly VectorToPoly(const DVector& kernelVector, const Basis<Mono>& startBasis){
	Poly ret;
	if(static_cast<size_t>(kernelVector.rows()) != startBasis.size()){
		std::cerr << "Error: the given Q column has " << kernelVector.rows()
			<< " rows, " << "but the given basis has " << startBasis.size() 
			<< " monomials. These must be the same." << std::endl;
		return ret;
	}
	for(auto row = 0; row < kernelVector.rows(); ++row){
		if(std::abs(kernelVector.coeff(row)) < EPSILON) continue;
		ret += kernelVector.coeff(row)*startBasis[row];
	}

	if(ret.size() == 0) return ret;
	return ret;
}

Poly ColumnToPoly(const Matrix& kernelMatrix, const Eigen::Index col, 
		const Basis<Mono>& startBasis){
	Poly ret;
	if(static_cast<size_t>(kernelMatrix.rows()) != startBasis.size()){
		std::cerr << "Error: the given Q matrix has " << kernelMatrix.rows()
			<< " rows, " << "but the given basis has " << startBasis.size() 
			<< " monomials. These must be the same." << std::endl;
		return ret;
	}
	for(Eigen::Index row = 0; row < kernelMatrix.rows(); ++row){
		if(kernelMatrix.coeff(row, col) == 0) continue;
		ret += kernelMatrix.coeff(row, col)*startBasis[row];
	}

	if(ret.size() == 0) return ret;
	coeff_class smallestCoeff = std::abs(ret[0].Coeff());
	for(auto& term : ret) smallestCoeff = std::min(std::abs(term.Coeff()), smallestCoeff);
	for(auto& term : ret) term /= smallestCoeff;
	return ret;
}

Poly ColumnToPoly(const DMatrix& kernelMatrix, const Eigen::Index col, 
		const Basis<Mono>& startBasis){
	Poly ret;
	if(static_cast<size_t>(kernelMatrix.rows()) != startBasis.size()){
		std::cerr << "Error: the given Q matrix has " << kernelMatrix.rows()
			<< " rows, " << "but the given basis has " << startBasis.size() 
			<< " monomials. These must be the same." << std::endl;
		return ret;
	}
	for(Eigen::Index row = 0; row < kernelMatrix.rows(); ++row){
		if(kernelMatrix.coeff(row, col) == 0) continue;
		ret += kernelMatrix.coeff(row, col)*startBasis[row];
	}

	if(ret.size() == 0) return ret;
	coeff_class smallestCoeff = std::abs(ret[0].Coeff());
	for(auto& term : ret) smallestCoeff = std::min(std::abs(term.Coeff()), smallestCoeff);
	for(auto& term : ret) term /= smallestCoeff;
	return ret;
}

DMatrix ExtractQMatrix(const Eigen::FullPivHouseholderQR<DMatrix>& solver, 
		               const int){
	return solver.matrixQ();
}

DMatrix ExtractQMatrix(const Eigen::ColPivHouseholderQR<DMatrix>& solver, 
		               const int dimension){
	return solver.householderQ()*DMatrix::Identity(dimension, dimension);
}

void ClearZeros(DMatrix* toClear) {
	if(!toClear){
		std::cerr << "Error: asked to clear the zeros from a nullptr instead of"
			<< " a matrix." << std::endl;
		return;
	}
	coeff_class threshold = EPSILON*toClear->cwiseAbs().maxCoeff();
	for(Eigen::Index row = 0; row < toClear->rows(); ++row){
		for(Eigen::Index col = 0; col < toClear->cols(); ++col){
			if(std::abs((*toClear)(row, col)) < threshold){
				(*toClear)(row, col) = 0;
			}
		}
	}
}

Basis<Mono> MinimalBasis(const std::vector<Poly>& polynomials) {
	Poly combinedPoly;
	for (const auto& poly : polynomials) {
		combinedPoly += poly;
	}
	std::vector<Mono> allUsedMonos(combinedPoly.size());
	for (auto i = 0u; i < combinedPoly.size(); ++i) {
		allUsedMonos[i] = combinedPoly[i]/combinedPoly[i].Coeff();
	}
	return Basis<Mono>(allUsedMonos);
}
