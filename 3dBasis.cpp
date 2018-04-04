#include "3dBasis.hpp"

int main(int argc, char* argv[]) {
    gsl_set_error_handler(&GSLErrorHandler);

    Arguments args = ParseArguments(argc, argv);
    if (args.options & OPT_VERSION){
        std::cout << "This is 3dBasis version " << VERSION << ", released "
            << RELEASE_DATE << ". The latest updates can always be found at "
            << "https://github.com/chussong/3dBasis." << std::endl;
        return EXIT_SUCCESS;
    }

    if (args.options & OPT_TEST) {
        return Test::RunAllTests();
    }

    if (args.degree == 0 || args.numP == 0){
        std::cerr << "Error: you must enter a number of particles, a degree, "
            << "and a value for delta." << std::endl;
        return EXIT_FAILURE;
    }

    // initialize all multinomials which might come up
    //
    // this is obviously something of a blunt instrument and could easily be
    // made more efficient
    for (int n = 1; n <= args.numP; ++n) {
        Multinomial::Initialize(n, 2*args.degree);
    }

    if (args.options & OPT_MULTINOMTEST) {
        for (char n = 0; n <= args.degree; ++n) {
            *args.outputStream << "n = " << std::to_string(n) << ": ";
            for (auto& mVector : Multinomial::GetMVectors(args.numP, n)) {
                //std::cout << std::endl << Multinomial::MVectorOut(mVector) << ": ";
                *args.outputStream << Multinomial::Lookup(args.numP, mVector) << ", ";
            }
            *args.outputStream << std::endl;
        }
    if (args.outputStream->rdbuf() != std::cout.rdbuf()) {
        delete args.outputStream;
    }
        return EXIT_SUCCESS;
    }

    if (args.options & OPT_STATESONLY) {
        ComputeBasisStates(args);
        if (args.outputStream->rdbuf() != std::cout.rdbuf()) {
            delete args.outputStream;
        }
        return EXIT_SUCCESS;
    }

    DMatrix hamiltonian = ComputeHamiltonian(args);
    if (args.outputStream->rdbuf() != std::cout.rdbuf()) {
        delete args.outputStream;
    }
    return EXIT_SUCCESS;
}

// return basis polynomials. They are NOT normalized w.r.t. partitions
std::vector<Poly> ComputeBasisStates(const Arguments& args) {
    int numP = args.numP;
    int degree = args.degree + args.numP; // add required Dirichlet derivatives
    int options = args.options;

    *args.outputStream << "(*Orthogonal basis states with N=" << numP << ", L="
            << degree << " (including Dirichlet derivatives).*)" << std::endl;
    
    std::vector<Basis<Mono>> allEvenBases;
    std::vector<Basis<Mono>> allOddBases;
    for (int deg = numP; deg <= degree; ++deg) {
            splitBasis<Mono> degBasis(numP, deg, options);
            allEvenBases.push_back(degBasis.EvenBasis());
            allOddBases.push_back(degBasis.OddBasis());
    }

    *args.outputStream << "(*EVEN STATE ORTHOGONALIZATION*)" << std::endl;
    std::vector<Poly> basisEven = ComputeBasisStates_SameParity(allEvenBases, args);

    *args.outputStream << "(*ODD STATE ORTHOGONALIZATION*)" << std::endl;
    std::vector<Poly> basisOdd = ComputeBasisStates_SameParity(allOddBases, args);

    *args.outputStream << std::endl;

    basisEven.insert(basisEven.end(), basisOdd.begin(), basisOdd.end());
    return basisEven;
}

// return basis polynomials. They are NOT normalized w.r.t. partitions
std::vector<Poly> ComputeBasisStates_SameParity(
        const std::vector<Basis<Mono>>& inputBases, const Arguments& args) {
    std::ostream& outStream = *args.outputStream;
    std::vector<Poly> orthogonalized = Orthogonalize(inputBases, outStream);

    // Basis<Mono> minimalBasis(MinimalBasis(orthogonalized));
    // if (outStream.rdbuf() == std::cout.rdbuf()) {
        // outStream << "Minimal basis: " << minimalBasis << std::endl;
    // } else {
        // outStream << "minimalBasis = " << MathematicaOutput(minimalBasis) 
            // << std::endl;
    // }

    return orthogonalized;
}

DMatrix PolysOnMinBasis(const Basis<Mono>& minimalBasis,
        const std::vector<Poly> orthogonalized, std::ostream&) {
    DMatrix polysOnMinBasis(minimalBasis.size(), orthogonalized.size());
    for (std::size_t i = 0; i < orthogonalized.size(); ++i) {
            polysOnMinBasis.col(i) = minimalBasis.DenseExpressPoly(
                            orthogonalized[i] );
    }

    // outStream << "(*Polynomials on this basis (as rows, not columns!):*)\n"
        // "polysOnMinBasis = " << MathematicaOutput(polysOnMinBasis.transpose()) 
        // << std::endl;

    return polysOnMinBasis;
}
    

DMatrix ComputeHamiltonian(const Arguments& args) {
    int numP = args.numP;
    int degree = args.degree + args.numP; // add required Dirichlet derivatives
    int options = args.options;

    *args.outputStream << "(*Matrix element test with N=" << numP << ", L="
            << degree << " (including Dirichlet derivatives).*)" << std::endl;
    
    std::vector<Basis<Mono>> allEvenBases;
    std::vector<Basis<Mono>> allOddBases;
    for(int deg = numP; deg <= degree; ++deg){
            splitBasis<Mono> degBasis(numP, deg, options);
            allEvenBases.push_back(degBasis.EvenBasis());
            allOddBases.push_back(degBasis.OddBasis());
    }

    *args.outputStream << "(*EVEN STATE ORTHOGONALIZATION*)" << std::endl;
    DMatrix evenHam = ComputeHamiltonian_SameParity(allEvenBases, args);

    *args.outputStream << "(*ODD STATE ORTHOGONALIZATION*)" << std::endl;
    DMatrix oddHam  = ComputeHamiltonian_SameParity(allOddBases, args);

    *args.outputStream << std::endl;

    return DMatrix();
}

// the portion of the hamiltonian computation which assumes equal pt parity
DMatrix ComputeHamiltonian_SameParity(const std::vector<Basis<Mono>>& inputBases,
                                      const Arguments& args) {
    Timer timer;
    std::ostream& outStream = *args.outputStream;
    std::vector<Poly> orthogonalized = ComputeBasisStates_SameParity(inputBases,
            args);

    Basis<Mono> minimalBasis(MinimalBasis(orthogonalized));
    if (outStream.rdbuf() == std::cout.rdbuf()) {
        outStream << "Minimal basis: " << minimalBasis << std::endl;
    } else {
        outStream << "minimalBasis = " << MathematicaOutput(minimalBasis) 
            << std::endl;
    }

    DMatrix polysOnMinBasis = PolysOnMinBasis(minimalBasis, orthogonalized,
            outStream);

    // outStream << "polysOnMinBasis:\n" << polysOnMinBasis << std::endl;
    DMatrix discPolys = DiscretizePolys(polysOnMinBasis, args.partitions);
    if (outStream.rdbuf() != std::cout.rdbuf()) {
        outStream << "(*Polynomials on this basis (as rows, not columns!):*)\n"
                << "polysOnMinBasis = " 
        << MathematicaOutput(polysOnMinBasis.transpose()) << std::endl;
        outStream << "(*And discretized:*)\ndiscretePolys = "
            << MathematicaOutput(discPolys.transpose()) << std::endl;
    }

    std::cout << "Fock space inner product for confirmation; monos:" << std::endl;
    DMatrix discGram(GramMatrix(minimalBasis, args.partitions, args.partWidth));
    DMatrix discGram_BasisStates = 
            discPolys.transpose() * discGram * discPolys;
    std::cout << discGram << std::endl << "basis states:" 
            << std::endl << discGram_BasisStates << std::endl;

    timer.Start();
    DMatrix monoMassMatrix(MassMatrix(minimalBasis, args.partitions,
                args.partWidth));
    if (outStream.rdbuf() != std::cout.rdbuf()) {
        outStream << "minBasisMassMatrix = "
            << MathematicaOutput(monoMassMatrix) << std::endl;
    } else {
        outStream << "Here is the monomial mass matrix we computed:\n"
            << monoMassMatrix << std::endl;
    }

    // DMatrix discMonoMass = PartitionMu_Mass(minimalBasis, monoMassMatrix, 
            // args.partitions, args.partWidth);
    // if (outStream.rdbuf() != std::cout.rdbuf()) {
        // outStream << "discretizedMinBasisMass = "
            // << MathematicaOutput(discMonoMass) << std::endl;
    // } else {
        // outStream << "Here it is discretized by mu:\n"
            // << discMonoMass << std::endl;
    // }

    DMatrix polyMassMatrix = discPolys.transpose()*monoMassMatrix*discPolys;
    if (outStream.rdbuf() != std::cout.rdbuf()) {
        outStream << "basisStateMassMatrix = "
                << MathematicaOutput(polyMassMatrix) << std::endl;
    } else {
        outStream << "Computed a mass matrix for the basis in " 
                << timer.TimeElapsedInWords() << ", getting this:\n"
                << polyMassMatrix << std::endl;
        EigenSolver solver(polyMassMatrix.cast<builtin_class>());
        outStream << "Here are the eigenvalues:\n" << solver.eigenvalues()
            << std::endl;
    }

    DMatrix monoKineticMatrix(KineticMatrix(minimalBasis, args.partitions,
                args.partWidth));
    if (outStream.rdbuf() != std::cout.rdbuf()) {
        outStream << "minBasisKineticMatrix = "
            << MathematicaOutput(monoKineticMatrix) << std::endl;
    } else {
        outStream << "Here is the monomial kinetic matrix we computed:\n"
            << monoKineticMatrix << std::endl;
    }

    DMatrix polyKineticMatrix = discPolys.transpose()*monoKineticMatrix*discPolys;
    if (outStream.rdbuf() != std::cout.rdbuf()) {
        outStream << "basisStateKineticMatrix = "
                << MathematicaOutput(polyKineticMatrix) << std::endl;
    } else {
        outStream << "Computed a kinetic matrix for the basis in " 
                << timer.TimeElapsedInWords() << ", getting this:\n"
                << polyKineticMatrix << std::endl;
        EigenSolver solver(polyKineticMatrix.cast<builtin_class>());
        outStream << "Here are the eigenvalues:\n" << solver.eigenvalues()
            << std::endl;
    }

    DMatrix hamiltonian = polyMassMatrix + polyKineticMatrix;
    if (outStream.rdbuf() != std::cout.rdbuf()) {
        outStream << "hamiltonian = "
                << MathematicaOutput(hamiltonian) << std::endl;
    } else {
        outStream << "Computed a Hamiltonian matrix for the basis in " 
                << timer.TimeElapsedInWords() << ", getting this:\n"
                << hamiltonian << std::endl;
        EigenSolver solver(hamiltonian.cast<builtin_class>());
        outStream << "Here are the eigenvalues:\n" << solver.eigenvalues()
            << std::endl;
    }


    return hamiltonian;
}

Arguments ParseArguments(int argc, char* argv[]) {
    std::vector<std::string> options;
    std::string arg;
    Arguments ret;
    ret.delta = 0;
    int j = 0;
    for(int i = 1; i < argc; ++i){
        arg = argv[i];
        if (arg.size() > 0) {
            if (arg[0] == '-') {
                if (arg.size() > 1 && arg[1] == 'o') {
                    // open next argument as outstream, appending to it
                    ret.outputStream = new std::ofstream(argv[i+1], 
                                    std::ios_base::out | std::ios_base::app);
                    ++i; // next argument is the filename so don't process it
                } else if (arg.size() > 1 && arg[1] == 'O') {
                    // open next argument as outstream, replacing it
                    ret.outputStream = new std::ofstream(argv[i+1], 
                                    std::ios_base::out | std::ios_base::trunc);
                    ++i; // next argument is the filename so don't process it
                } else {
                    options.push_back(arg);
                }
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
                            std::cerr << "Error: at most three non-option arguments"
                                    << " may be given." << std::endl;
                            return ret;
                }
                ++j;
            }
        }
    }
    if(j < 2) ret.numP = 0; // invalidate the input since it was insufficient
    ret.options = ParseOptions(options);
    if (argc < 3 || std::abs<builtin_class>(ret.delta) < EPSILON) ret.delta = 0.5;
    return ret;
}

int ParseOptions(std::vector<std::string> options) {
    int ret = 0;
    for(auto& opt : options){
        if(opt.compare(0, 2, "-d") == 0){
            ret |= OPT_DEBUG;
            ret |= OPT_OUTPUT;
            continue;
        }
        if(opt.compare(0, 2, "-i") == 0){
            ret |= OPT_IPTEST;
            continue;
        }
        if(opt.compare(0, 2, "-m") == 0){
            ret |= OPT_MULTINOMTEST;
            continue;
        }
        if(opt.compare(0, 2, "-M") == 0){
            ret |= OPT_ALLMINUS;
            continue;
        }
        if(opt.compare(0, 2, "-o") == 0){
            ret |= OPT_OUTPUT;
            continue;
        }
        if(opt.compare(0, 2, "-s") == 0){
            ret |= OPT_STATESONLY;
            continue;
        }
        if(opt.compare(0, 2, "-t") == 0){
            ret |= OPT_TEST;
            continue;
        }
        if(opt.compare(0, 2, "-v") == 0){
            ret |= OPT_VERSION;
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

void GSLErrorHandler(const char* reason, const char* file, int line, int err) {
    std::cerr << "GSL Error in " << file << ":" << line << " --- "
        << gsl_strerror(err) << ", " << reason << std::endl;
}

bool particle::operator==(const particle& other) const {
	return (pm == other.pm) && (pt == other.pt);
}

/*
std::list<Triplet> ConvertToRows(const std::vector<Poly>& polyForms, 
		const Basis<Mono>& targetBasis, const Eigen::Index rowOffset) {
	if(polyForms.size() == 0) return std::list<Triplet>();
	std::list<Triplet> ret = targetBasis.ExpressPoly(polyForms[0], 0,
			rowOffset);
	for(auto i = 1u; i < polyForms.size(); ++i){
		ret.splice(ret.end(), targetBasis.ExpressPoly(polyForms[i], i,
				rowOffset));
	}
	return ret;
}
*/

/*Poly ColumnToPoly(const SMatrix& kernelMatrix, const Eigen::Index col, 
		const Basis<Mono>& startBasis) {
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
}*/

Poly ColumnToPoly(const DMatrix& kernelMatrix, const Eigen::Index col, 
		const Basis<Mono>& startBasis) {
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
    // coeff_class smallestCoeff = std::abs(ret[0].Coeff());
    // for(auto& term : ret) smallestCoeff = std::min(std::abs(term.Coeff()), smallestCoeff);
    // for(auto& term : ret) term /= smallestCoeff;
    return ret;
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
            if(std::abs<builtin_class>((*toClear)(row, col)) < threshold){
                    (*toClear)(row, col) = 0;
            }
        }
    }
}
