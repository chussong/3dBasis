#include "3dBasis.hpp"

int main(int argc, char* argv[]) {
	arguments args = ParseArguments(argc, argv);
	if(args.options & OPT_VERSION){
		std::cout << "This is 3dBasis version " << VERSION << ", released "
			<< RELEASE_DATE << ". The latest updates can always be found at "
			<< "https://github.com/chussong/3dBasis." << std::endl;
		return EXIT_SUCCESS;
	}

	if(args.degree == 0 || args.numP == 0){
		std::cerr << "Error: you must enter a number of particles, a degree, "
			<< "and a value for delta." << std::endl;
		return EXIT_FAILURE;
	}

	if(std::abs(args.delta) < EPSILON) args.delta = 0.5;

	if(args.options & OPT_IPTEST){
		return InnerProductTest(args);
	}

	if(args.options & OPT_MSORTING){
		args.options = args.options | OPT_EQNMOTION;
		FindPrimariesByM(args);
		// basis = Orthogonalize(primaries);
	}

	if(args.options & OPT_BRUTE){
		return FindPrimariesBruteForce(args);
	} else {
		return FindPrimariesParityOnly(args);
	}
}

//------------------------------------------------------------------------------
// targetBasis is not split because this would require finding the overlap
// between two kernels with distinct bases, which is quite nontrivial. It's 
// possible that splitting targetBasis could fit well with a separate 
// orthogonalization step if one is added in the future.
//------------------------------------------------------------------------------
int FindPrimaries(const arguments& args){
	int numP = args.numP;
	int degree = args.degree;
	coeff_class delta = args.delta;
	int options = args.options;
	splitBasis<mono> startingBasisA(numP, degree, options);
	splitBasis<mono> startingBasisS = startingBasisA.BecomeAsymmetric();
	std::cout << "Constructed an asymmetric starting basis with " 
		<< startingBasisA.EvenBasis().size()
		<< " even elements and " << startingBasisA.OddBasis().size()
		<< " odd elements." << std::endl;
	std::cout << "Constructed a symmetric starting basis with " 
		<< startingBasisS.EvenBasis().size()
		<< " even elements and " << startingBasisS.OddBasis().size()
		<< " odd elements." << std::endl;
	Basis<mono> targetBasis(numP, degree-1, options);
	std::cout << "Constructed a target basis with " << targetBasis.size()
		<< " elements." << std::endl;

	// - create matrix of K acting on each element of startingBasis
	Matrix evenKActionA(KMatrix(startingBasisA.EvenBasis(), targetBasis, delta,
				options));
	Matrix oddKActionA (KMatrix(startingBasisA.OddBasis() , targetBasis, delta,
				options));
	Matrix evenKActionS(KMatrix(startingBasisS.EvenBasis(), targetBasis, delta,
				options));
	Matrix oddKActionS (KMatrix(startingBasisS.OddBasis() , targetBasis, delta,
				options));

	// - find kernel of above matrix and output

	std::vector<poly> evenKernelA = Kernel(evenKActionA, startingBasisA.EvenBasis(),
			options, false);
	std::vector<poly> oddKernelA  = Kernel(oddKActionA , startingBasisA.OddBasis(),
			options, false);
	std::vector<poly> evenKernelS = Kernel(evenKActionS, startingBasisS.EvenBasis(),
			options, false);
	std::vector<poly> oddKernelS  = Kernel(oddKActionS , startingBasisS.OddBasis(),
			options, false);

	std::cout << "Found a total of " 
		<< 2*evenKernelA.size() + 2*oddKernelA.size()\
			+ evenKernelS.size() + oddKernelS.size();
	if(options & OPT_DIRICHLET) std::cout << " Dirichlet";
	std::cout << " primaries." << std::endl;

	/*std::cout << "Even asymmetric:" << std::endl;
	for(auto& kernelVec : evenKernelA) std::cout << kernelVec << std::endl;
	std::cout << "Odd asymmetric:" << std::endl;
	for(auto& kernelVec : oddKernelA) std::cout << kernelVec << std::endl;
	std::cout << "Even symmetric:" << std::endl;
	for(auto& kernelVec : evenKernelS) std::cout << kernelVec << std::endl;
	std::cout << "Odd symmetric:" << std::endl;
	for(auto& kernelVec : oddKernelS) std::cout << kernelVec << std::endl;*/
	
	return EXIT_SUCCESS;
}

int FindPrimariesParityOnly(const arguments& args){
	int numP = args.numP;
	int degree = args.degree;
	coeff_class delta = args.delta;
	int options = args.options;

	//options = options | OPT_DEBUG;

	splitBasis<mono> startingBasis(numP, degree, options);
	std::cout << "Constructed a starting basis with " 
		<< startingBasis.EvenBasis().size()
		<< " even elements and " << startingBasis.OddBasis().size()
		<< " odd elements." << std::endl;
	if(options & OPT_OUTPUT){
		std::cout << startingBasis.EvenBasis() << std::endl;
		std::cout << startingBasis.OddBasis() << std::endl;
	}
	//WARNING: BELOW LINE IS TEMPORARY ONLY
	options = options & (~OPT_DIRICHLET);
	Basis<mono> targetBasis(numP, degree-1, options);
	std::cout << "Constructed a target basis with " << targetBasis.size()
		<< " elements." << std::endl;
	if(options & OPT_OUTPUT){
		std::cout << targetBasis << std::endl;
	}

	Matrix evenKAction(KMatrix(startingBasis.EvenBasis(), targetBasis, delta,
				options));
	Matrix oddKAction (KMatrix(startingBasis.OddBasis() , targetBasis, delta,
				options));

	std::vector<poly> evenKernel = Kernel(evenKAction, startingBasis.EvenBasis(),
			options, options & OPT_OUTPUT);
	std::vector<poly> oddKernel  = Kernel(oddKAction , startingBasis.OddBasis(),
			options, options & OPT_OUTPUT);

	std::cout << "Found a total of " << evenKernel.size() << " even and "
		<< oddKernel.size() << " odd ";
	if(options & OPT_DIRICHLET) std::cout << " Dirichlet";
	std::cout << "primaries." << std::endl;
	if(options & OPT_OUTPUT){
		std::cout << "Even:" << std::endl;
		std::cout << evenKernel << std::endl;
		/*for(auto& kernelVec : evenKernel){
			std::cout << kernelVec.HumanReadable() << std::endl;
		}*/
		std::cout << "Odd:" << std::endl;
		std::cout << oddKernel << std::endl;
		/*for(auto& kernelVec : oddKernel){
			std::cout << kernelVec.HumanReadable() << std::endl;
		}*/
	}

	return EXIT_SUCCESS;
}

int FindPrimariesBruteForce(const arguments& args){
	int numP = args.numP;
	int degree = args.degree;
	coeff_class delta = args.delta;
	int options = args.options;

	//options = options | OPT_DEBUG;

	Basis<mono> startingBasis(numP, degree, options);
	std::cout << "Constructed a starting basis with " << startingBasis.size()
		<< " elements." << std::endl;
	if(options & OPT_DEBUG){
		std::cout << startingBasis << std::endl;
	}
	Basis<mono> targetBasis(numP, degree-1, options);
	std::cout << "Constructed a target basis with " << targetBasis.size()
		<< " elements." << std::endl;
	if(options & OPT_DEBUG){
		std::cout << targetBasis << std::endl;
	}

	// - create matrix of K acting on each element of startingBasis
	Matrix kAction(KMatrix(startingBasis, targetBasis, delta, options));

	// - find kernel of above matrix and output

	std::vector<poly> kernel = Kernel(kAction, startingBasis, options,
			options & OPT_DEBUG);

	if(options & OPT_DEBUG){
		std::cout << "Found the following " << kernel.size() << "-dimensional kernel:"
			<< std::endl;
		for(auto& kernelVec : kernel) std::cout << kernelVec.HumanReadable() << std::endl;
	} else {
		std::cout << "Found " << kernel.size();
		if(options & OPT_DIRICHLET) std::cout << " Dirichlet";
		std::cout << " primary operators." << std::endl;
	}
	
	return EXIT_SUCCESS;
}

int FindPrimariesByM(const arguments& args){
	mBasis startBasis(args.numP, args.degree, args.options);
	mBasis targetBasis(args.numP, args.degree-1, args.options);

	std::vector<poly> primaries;
	unsigned int primCount = 0u;
	// I think this is supposed to be -= 2, but maybe -= 1 is needed
	for(int L = args.degree; L >= 0; L -= 1){
		primCount += AddPrimariesAtL(startBasis, targetBasis, L, primaries,
				args.options, args.delta);
	}

	std::cout << primCount << " primaries identified: " << primaries
		<< std::endl;

	return EXIT_SUCCESS;
}

unsigned int AddPrimariesAtL(const mBasis& startBasis, const mBasis& targetBasis,
		const unsigned int L, std::vector<poly>& primaries, 
		const coeff_class delta, const int options){
	//std::cout << "Constructing polynomials of L=" << L << " actions..." << std::endl;
	std::vector<poly> K1Actions, K2Actions, K3Actions;
	for(auto& topState : startBasis.TopStates(L)){
		if(L+1 < targetBasis.Degree()) K1Actions.emplace_back(topState.K1(delta));
		if(L   < targetBasis.Degree()) K2Actions.emplace_back(topState.K2(delta));
		K3Actions.emplace_back(topState.K3(delta));
	}

	//std::cout << "Converting L actions to triplets..." << std::endl;
	std::list<Triplet> entries;
	if(L+1 <= targetBasis.Degree()){
		/*std::cout << "Querying targetBasis.Level(" << L+1 << "); the basis has "
			<< "degree " << targetBasis.Degree() << ", and the size of the "
			<< "level in question is " << targetBasis.LevelSize(L+1) << "." 
			<< std::endl;
		std::cout << *targetBasis.Level(L+1) << std::endl;
		std::cout << "Meanwhile, the matrix we're using is this:" << std::endl;
		std::cout << K1Actions << std::endl;*/
		entries.splice(entries.end(), 
				ConvertToRows(K1Actions, *targetBasis.Level(L+1), 0));
	}
	if(L <= targetBasis.Degree()){
		entries.splice(entries.end(), 
				ConvertToRows(K2Actions, *targetBasis.Level(L), 
					targetBasis.LevelSize(L+1)));
	}
	entries.splice(entries.end(), 
			ConvertToRows(K3Actions, *targetBasis.Level(L-1), 
				targetBasis.LevelSize(L+1) + targetBasis.LevelSize(L)));

	//std::cout << "List of triplets done, matrixifying..." << std::endl;
	Matrix matrixK(targetBasis.LevelSize(L+1) + targetBasis.LevelSize(L) + 
			targetBasis.LevelSize(L-1), startBasis.TopStates(L).size());

	// if assignment is slow, can do ret.reserve(3*numP) to speed it up
	matrixK.setFromTriplets(entries.begin(), entries.end());
	// matrix must be compressed here but setFromTriplets does it automatically

	std::cout << "Computing kernel of K matrix at level " << L << "." << std::endl;
	std::vector<poly> kernel = Kernel(matrixK, 
			Basis<poly>(startBasis.TopStates(L)), options, false);
	for(auto& newPrimary : kernel) primaries.push_back(std::move(newPrimary));
	// if the pushing back is slow, can easily resize and fill in
	//std::cout << "Kernel:" << std::endl;
	//std::cout << kernel << std::endl;
	return kernel.size()*(2*L+1);
}

int InnerProductTest(const arguments& args){
	int numP = args.numP;
	int degree = args.degree;
	//coeff_class delta = args.delta;
	int options = args.options;

	//options = options | OPT_DEBUG;

	std::cout << "Beginning inner product test with N=" << numP << ", L="
		<< degree << "." << std::endl;
	splitBasis<mono> basis(numP, degree, options);
	/*std::cout << "Even with even: " << std::endl;
	for(auto& monoA : basis.EvenBasis()){
		for(auto& monoB : basis.EvenBasis()){
			std::cout << monoA << " x " << monoB << " = " 
				<< mono::InnerProduct(monoA, monoB) << std::endl;
		}
	}*/
	/*std::cout << "Even with odd: " << std::endl;
	for(auto& monoA : basis.EvenBasis()){
		for(auto& monoB : basis.OddBasis()){
			std::cout << monoA << " x " << monoB << " = " 
				<< mono::InnerProduct(monoA, monoB) << std::endl;
		}
	}
	std::cout << "Odd with even: " << std::endl;
	for(auto& monoA : basis.OddBasis()){
		for(auto& monoB : basis.EvenBasis()){
			std::cout << monoA << " x " << monoB << " = " 
				<< mono::InnerProduct(monoA, monoB) << std::endl;
		}
	}*/
	/*std::cout << "Odd with odd: " << std::endl;
	for(auto& monoA : basis.OddBasis()){
		for(auto& monoB : basis.OddBasis()){
			std::cout << monoA << " x " << monoB << " = " 
				<< mono::InnerProduct(monoA, monoB) << std::endl;
		}
	}*/

	std::vector<Basis<mono>> allEvenBases;
	std::vector<Basis<mono>> allOddBases;
	for(int deg = numP; deg <= degree; ++deg){
		splitBasis<mono> degBasis(numP, deg, options);
		allEvenBases.push_back(degBasis.EvenBasis());
		allOddBases.push_back(degBasis.OddBasis());
	}

	std::cout << "EVEN STATE ORTHOGONALIZATION" << std::endl;
	Orthogonalize(allEvenBases);

	std::cout << "ODD STATE ORTHOGONALIZATION" << std::endl;
	Orthogonalize(allOddBases);

	return EXIT_SUCCESS;
}

int Orthogonalize(const std::vector<Basis<mono>>& inputBases){
	Timer timer;
	DMatrix gram = GramMatrix(inputBases);
	std::cout << "Gram matrix constructed in " << timer.TimeElapsedInWords()
		<< "." << std::endl;
	if(TotalSize(inputBases) <= 7){
		std::cout << gram << std::endl;
	}

	timer.Start();
	DQRSolver solver;
	// it's possible there's a smarter way to set the threshold than this,
	// which could possibly take a long time? We multiply our custom epsilon by
	// the largest coefficient of gram's component-wise absolute value
	solver.setThreshold(EPSILON*gram.cwiseAbs().maxCoeff());
	solver.compute(gram);
	DMatrix QMatrix = ExtractQMatrix(solver, gram.rows());
	std::cout << "Q matrix found in " << timer.TimeElapsedInWords() << "."
		<< std::endl;
	if(TotalSize(inputBases) <= 7){
		std::cout << QMatrix << std::endl;
	}
	std::cout << "Rank, i.e. number of independent operators: " << solver.rank() 
		<< std::endl;
	return solver.rank();
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
		if(opt.compare(0, 2, "-b") == 0){
			ret = ret | OPT_BRUTE;
			continue;
		}
		if(opt.compare(0, 2, "-p") == 0){
			ret = ret | OPT_PARITYONLY;	// currently the default behavior
			continue;
		}
		if(opt.compare(0, 2, "-e") == 0){
			ret = ret | OPT_EQNMOTION;
			continue;
		}
		if(opt.compare(0, 2, "-m") == 0){
			ret = ret | OPT_MSORTING;
			continue;
		}
		if(opt.compare(0, 2, "-d") == 0){
			ret = ret | OPT_DEBUG;
			ret = ret | OPT_OUTPUT;
			continue;
		}
		if(opt.compare(0, 2, "-D") == 0){
			ret = ret | OPT_DIRICHLET;
			ret = ret | OPT_EQNMOTION;
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
	}
	return ret;
}

bool particle::operator==(const particle& other) const{
	return (pm == other.pm) && (pt == other.pt) && (pp == other.pp);
}

std::vector<Basis<mono>> mBasis::MakeNegativeLevels(
		std::vector<Basis<mono>>& nonNegativeMLevels){
	std::vector<Basis<mono>> ret;
	std::vector<mono> newLevel;
	for(auto degree = nonNegativeMLevels.size()-1; degree > 0; --degree){
		newLevel.clear();
		for(auto& mon : nonNegativeMLevels[degree]){
			newLevel.push_back(mon.MirrorPM());
		}
		ret.emplace_back(newLevel);
	}
	for(auto degree = 0u; degree < nonNegativeMLevels.size(); ++degree){
		ret.push_back(nonNegativeMLevels[degree]);
	}
	return ret;
}

mBasis::mBasis(const int numP, const int degree, const int options){
	std::vector<Basis<mono>> nonNegativeMLevels;
	for(int M = 0; M <= degree; ++M){
		nonNegativeMLevels.push_back(BasisAtM(numP, degree, M, options));
		if(options & OPT_DEBUG) std::cout << nonNegativeMLevels[M] << std::endl;
	}
	mLevels = MakeNegativeLevels(nonNegativeMLevels);

	Matrix L3Actions;
	std::vector<poly> L3Kernel;
	for(int M = 0; M < degree; ++M){
		L3Actions = L3Matrix(*Level(M), *Level(M+1));
		L3Kernel = Kernel(L3Actions, *Level(M), options, true);
		//std::cout << "States at the top of the L=" << M << " multiplet: " 
			//<< L3Kernel << std::endl;

		topStates.push_back(L3Kernel);
	}
	// all states with M=d are in L3Kernel;
	L3Kernel.resize(LevelSize(degree));
	for(auto i = 0u; i < LevelSize(degree); ++i){
		L3Kernel[i] = poly((*Level(degree))[i]);
	}
	topStates.push_back(L3Kernel);
}

// This function was written before I figured out that the EoM were necessary
// for an L eigenbasis, so it generates non-EoM compliant states and then throws
// them out. If you care more than I do right now, you can change it to directly
// generate only the useful states.
Basis<mono> mBasis::BasisAtM(const int numP, const int degree, const int M, 
		const int options){
	if(options & OPT_DEBUG) 
		std::cout << "Constructing basis at M = " << M << "." << std::endl;
	std::vector<std::vector<int>> minusCfgs = GetStatesUpToDegree(numP,
			(degree + M)/2, M);
	int zeroCount, totalPm;
	std::vector<std::vector<int>> plusCfgs, perpCfgs;
	std::vector<std::vector<particle>> combinedCfgs;
	std::vector<particle> combinedCfg;
	combinedCfg.resize(numP);
	for(auto& minusCfg : minusCfgs){
		zeroCount = numP;
		totalPm = 0;
		plusCfgs.clear();
		while(zeroCount > 0){
			if(minusCfg[numP-zeroCount] != 0){
				totalPm += minusCfg[numP-zeroCount];
				--zeroCount;
			} else {
				break;
			}
		}
		if(zeroCount == 0 && totalPm > M) continue;
		plusCfgs = GetStatesAtDegree(zeroCount, totalPm - M);
		for(auto& plusCfg : plusCfgs){
			plusCfg.insert(plusCfg.begin(), numP - plusCfg.size(), 0);
			if(options & OPT_DEBUG) std::cout << "Minus: " << minusCfg << std::endl;
			if(options & OPT_DEBUG) std::cout << "Plus: " << plusCfg << std::endl;
			//combinedCfg.clear();
			for(auto i = 0; i < numP; ++i){
				combinedCfg[i].pm = minusCfg[i];
				combinedCfg[i].pp = plusCfg[i];
			}
			combinedCfgs.push_back(combinedCfg);
		}
	}

	std::vector<mono> basisMonos;
	std::vector<int> nodes;
	//const bool useEoM = (options & OPT_EQNMOTION);
	const bool useEoM = true; // EoM is mandatory for L eigenbasis
	for(auto& cfg : combinedCfgs){
		if(options & OPT_DEBUG) std::cout << "Processing this cfg: " << cfg << std::endl;
		nodes = IdentifyNodes(cfg);
		int remainingEnergy = degree;
		for(auto& part : cfg) remainingEnergy -= part.pm + part.pp;
		if(options & OPT_DEBUG){
			std::cout << "Generating states with remaining energy " << remainingEnergy
			<< std::endl;
		}
		std::vector<std::vector<int>> perp(CfgsFromNodes(remainingEnergy, nodes,
															true));
		for(auto& newCfg : CombinedCfgs(cfg, perp, 2)){
			basisMonos.emplace_back(newCfg, useEoM);
		}
	}

	/*for(mono& m : basisMonos){
		m.Order();
		m.Coeff() = 1;
	}*/

	return Basis<mono>(basisMonos);
}

Matrix mBasis::L3Matrix(const Basis<mono>& startingMBasis, const Basis<mono>& targetMBasis){
	//std::cout << "Constructing polynomials of L actions..." << std::endl;
	if(startingMBasis.size() == 0) return Matrix();
	std::vector<poly> L3Actions;
	for(auto& basisMono : startingMBasis){
		L3Actions.emplace_back(basisMono.L3());
	}

	//std::cout << "Converting L actions to triplets..." << std::endl;
	std::list<Triplet> triplets = ConvertToRows(L3Actions, targetMBasis, 0);

	//std::cout << "List of triplets done, matrixifying..." << std::endl;
	Matrix ret(targetMBasis.size(), startingMBasis.size());
	// if this is slow, can do ret.reserve(3*numP) to speed it up
	ret.setFromTriplets(triplets.begin(), triplets.end());
	// matrix must be compressed here but setFromTriplets does it automatically

	//std::cout << "Matrix done, here it is:" << std::endl << ret << std::endl;
	return ret;
}

/*std::vector<poly> mBasis::CompleteMultiplet(const poly& topState){
	std::vector<poly> ret;
	ret.push_back(topState);
	poly nextState = topState.L1();
	while(nextState.size() != 0){
		ret.push_back(nextState);
		nextState = nextState.L1();
	}
	return ret;
}*/

const Basis<mono>* mBasis::Level(const int M) const {
	if(static_cast<unsigned int>(std::abs(M)) > Degree()){
		std::cerr << "Error: basis at M=" << M << " was requested, but only "
			<< mLevels.size() << " are known." << std::endl;
		return nullptr;
	}
	return &mLevels[M + Degree()];
}

unsigned int mBasis::LevelSize(const int M) const {
	if(static_cast<unsigned int>(std::abs(M)) > Degree()) return 0;
	return Level(M)->size();
}

// note: triplets displayed (row, column, value) despite matrix's storage type
std::ostream& operator<<(std::ostream& os, const Triplet& out){
	return os << "(" << out.row() << "," << out.col() << "," << out.value()
		<< ")";
}

/*
// note: I recommend starting a new line before this, but I'm not your parents
std::ostream& operator<<(std::ostream& os, const Matrix& out){
	out.isVector() ? os << "( " : os << "/ ";
	for(auto col = 0; col < out.cols(); ++col){
		os << out.coeff(0, col) << ", ";
	}
	out.isVector() ? return os << "\b\b )" : os << "\b\b \\\n";

	os << "| ";
	for(auto row = 1; row < out.rows()-1; ++row){
		for(auto col = 0; col < out.cols(); ++col){
			os << out.coeff(row, col) << ", ";
		}
	}
	os << "\b\b |\n";

	os << "\\ ";
	for(auto col = 0; col < out.cols(); ++col){
		of << out.coeff(out.rows()-1, col) << ", ";
	}
	return os << "\b\b /";
}*/

// this version is a 'loose' EoM compliance that only removes Pp which are on
// the same particle as a Pm
bool EoMAllowed(const std::vector<particle>& cfg){
	for(auto& p : cfg){
		if(p.pm != 0 && p.pp != 0) return false;
	}
	return true;
}

// this is a 'strict' EoM compliance where all Pp must be 0
/*bool EoMAllowed(const std::vector<particle>& cfg){
	for(auto& p : cfg){
		if(p.pp != 0) return false;
	}
	return true;
}*/

Matrix KMatrix(const Basis<mono>& startingBasis, const Basis<mono>& targetBasis,
		const coeff_class delta, const int options){
	//std::cout << "Constructing polynomials of K actions..." << std::endl;
	if(startingBasis.size() == 0) return Matrix(0, 0);
	std::vector<poly> K1Actions, K2Actions, K3Actions;
	for(auto& basisMono : startingBasis){
		K1Actions.emplace_back(basisMono.K1(delta));
		K2Actions.emplace_back(basisMono.K2(delta));
		if(options & OPT_DIRICHLET){
			K3Actions.push_back(poly::DeleteNonDirichlet(basisMono.K3(delta)));
		} else {
			K3Actions.emplace_back(basisMono.K3(delta));
		}
	}

	//std::cout << "Converting K actions to triplets..." << std::endl;
	std::list<Triplet> entries = ConvertToRows(K1Actions, targetBasis, 0);
	entries.splice(entries.end(), ConvertToRows(K2Actions, targetBasis, 
				targetBasis.size()));
	entries.splice(entries.end(), ConvertToRows(K3Actions, targetBasis, 
				2*targetBasis.size()));

	//std::cout << "List of triplets done, matrixifying..." << std::endl;
	Matrix ret(3*targetBasis.size(), startingBasis.size());
	// if this is slow, can do ret.reserve(3*numP) to speed it up
	ret.setFromTriplets(entries.begin(), entries.end());
	// matrix must be compressed here but setFromTriplets does it automatically
	return ret;
}

Matrix K13Matrix(const Basis<mono>& startingBasis, const Basis<mono>& targetBasis,
		const coeff_class delta){
	std::vector<poly> K1Actions, K3Actions;
	for(auto& basisMono : startingBasis){
		K1Actions.emplace_back(basisMono.K1(delta));
		K3Actions.emplace_back(basisMono.K3(delta));
	}

	std::list<Triplet> entries = ConvertToRows(K1Actions, targetBasis, 0);
	entries.splice(entries.end(), ConvertToRows(K3Actions, targetBasis, 
				targetBasis.size()));

	Matrix ret(2*targetBasis.size(), startingBasis.size());
	ret.setFromTriplets(entries.begin(), entries.end());
	return ret;
}

Matrix K2Matrix(const Basis<mono>& startingBasis, const Basis<mono>& targetBasis,
		const coeff_class delta){
	std::vector<poly> K2Actions;
	for(auto& basisMono : startingBasis){
		K2Actions.emplace_back(basisMono.K2(delta));
	}

	std::list<Triplet> entries = ConvertToRows(K2Actions, targetBasis, 0);

	Matrix ret(targetBasis.size(), startingBasis.size());
	ret.setFromTriplets(entries.begin(), entries.end());
	return ret;
}

std::array<Matrix,4> KMatrices(const splitBasis<mono>& startingBasis,
		const splitBasis<mono>& targetBasis, const coeff_class delta){
	return std::array<Matrix,4>({{
			K13Matrix(startingBasis.EvenBasis(), targetBasis.EvenBasis(), delta),
			K2Matrix (startingBasis.EvenBasis(), targetBasis.OddBasis(), delta),
			K13Matrix(startingBasis.OddBasis() , targetBasis.OddBasis(), delta),
			K2Matrix (startingBasis.OddBasis() , targetBasis.EvenBasis(), delta)
			}});
}

std::list<Triplet> ConvertToRows(const std::vector<poly>& polyForms, 
		const Basis<mono>& targetBasis, const Eigen::Index rowOffset){
	if(polyForms.size() == 0) return std::list<Triplet>();
	std::list<Triplet> ret = targetBasis.ExpressPoly(polyForms[0], 0,
			rowOffset);
	for(auto i = 1u; i < polyForms.size(); ++i){
		ret.splice(ret.end(), targetBasis.ExpressPoly(polyForms[i], i,
				rowOffset));
	}
	return ret;
}

// takes QR decomposition of the matrix and returns the polynomial forms of its
// rightmost N columns, which are the N orthonormal basis vectors of the kernel
/*std::vector<poly> Kernel(const Matrix& KActions, const Basis<mono>& startBasis,
		const bool outputKernel){
	if(KActions.rows() == 0 || KActions.cols() == 0) return std::vector<poly>();
	std::cout << "Computing kernel from K matrix..." << std::endl;
	//std::cout << KActions << std::endl;
	QRSolver solver;
	solver.setPivotThreshold(EPSILON); // norms smaller than this are zero
	solver.compute(KActions.transpose());
	std::cout << "Solved. Found rank " << solver.rank() << ", i.e. "
		<< startBasis.size() - solver.rank() << " kernel elements." << std::endl;

	std::cout << "Converting the kernel to polynomials..." << std::endl;


	DVector projector = Eigen::VectorXd::Zero(startBasis.size());
	DVector kernelVector(startBasis.size());
	std::vector<poly> ret;
	ret.resize(startBasis.size() - solver.rank());

	if(!outputKernel) return ret;

	for(auto col = 0u; col < startBasis.size() - solver.rank(); ++col){
		projector(solver.rank() + col-1) = 0;
		projector(solver.rank() + col) = 1;
		kernelVector = solver.matrixQ()*projector;
		//std::cout << "Projecting out with this: " << projector << std::endl;
		//std::cout << kernelVector << "\n----------" << std::endl;
		ret[col] = VectorToPoly(kernelVector, startBasis);
	}

	return ret;
}*/

std::vector<poly> CombineKernels(const std::vector<poly>& kernel1,
		const std::vector<poly>& kernel2){
	std::vector<poly> ret;
	for(auto& poly1 : kernel1){
		for(auto& poly2 : kernel2){
			if(poly1 == poly2){
				ret.push_back(poly1);
				break;
			}
		}
	}
	std::cout << "Combined kernels to get " << ret.size() << " primaries." << std::endl;
	return ret;
}

/*poly VectorToPoly(const Vector& kernelVector, const Basis<mono>& startBasis){
	poly ret;
	if(static_cast<size_t>(kernelMatrix.rows()) != startBasis.size()){
		std::cerr << "Error: the given Q matrix has " << kernelMatrix.rows()
			<< " rows, " << "but the given basis has " << startBasis.size() 
			<< " monomials. These must be the same." << std::endl;
		return ret;
	}
	for(auto row = 0; row < kernelVector.rows(); ++row){
		if(kernelVector.coeff(row, col) == 0) continue;
		ret += kernelMatrix.coeff(row, col)*startBasis[row];
	}

	if(ret.size() == 0) return ret;
	coeff_class smallestCoeff = std::abs(ret[0].Coeff());
	for(auto& term : ret) smallestCoeff = std::min(std::abs(term.Coeff()), smallestCoeff);
	for(auto& term : ret) term /= smallestCoeff;
	return ret;
}*/

// NOTE: THIS IS BROKEN
// We need to either find a way to put a [] operator (or equivalent) into the
// generic basis or we need to specialize this for the mono and poly variants
poly VectorToPoly(const DVector& kernelVector, const Basis<mono>& startBasis){
	poly ret;
	if(static_cast<size_t>(kernelVector.rows()) != startBasis.size()){
		std::cerr << "Error: the given Q column has " << kernelVector.rows()
			<< " rows, " << "but the given basis has " << startBasis.size() 
			<< " monomials. These must be the same." << std::endl;
		return ret;
	}
	for(auto row = 0; row < kernelVector.rows(); ++row){
		if(std::abs(kernelVector.coeff(row)) < 1e-10) continue;
		ret += kernelVector.coeff(row)*startBasis[row];
	}

	if(ret.size() == 0) return ret;
	coeff_class smallestCoeff = std::abs(ret[0].Coeff());
	for(auto& term : ret) smallestCoeff = std::min(std::abs(term.Coeff()), smallestCoeff);
	for(auto& term : ret) term /= smallestCoeff;
	return ret;
}

poly ColumnToPoly(const Matrix& kernelMatrix, const Eigen::Index col, 
		const Basis<mono>& startBasis){
	poly ret;
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

poly ColumnToPoly(const DMatrix& kernelMatrix, const Eigen::Index col, 
		const Basis<mono>& startBasis){
	poly ret;
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

/*DMatrix GramMatrix(const Basis<mono>& basis){
	std::vector<Triplet> entries;
	for(auto row = 0u; row < basis.size(); ++row){
		for(auto col = row; col < basis.size(); ++col){
			entries.emplace_back(row, col, 
					mono::InnerProduct(basis[row], basis[col]));
			if(row != col){
				entries.emplace_back(col, row, 
						mono::InnerProduct(basis[row], basis[col]));
			}
		}
	}
	Matrix gram(basis.size(), basis.size());
	gram.setFromTriplets(entries.begin(), entries.end());
	return gram;
}*/
