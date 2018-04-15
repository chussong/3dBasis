#include "calculation.hpp"

int Calculate(const Arguments& args) {
    gsl_set_error_handler(&GSLErrorHandler);

    if (args.options & OPT_TEST) {
        return Test::RunAllTests();
    }

    // initialize all multinomials which might come up
    //
    // this is obviously something of a blunt instrument and could easily be
    // made more efficient

    for (int n = 1; n <= args.numP; ++n) {
        std::cout << "Initialize(" << n << ", " << 2*args.degree << ")" 
            << std::endl;
        Multinomial::Initialize(n, 2*args.degree);
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

    if (outStream.rdbuf() != std::cout.rdbuf()) {
        DMatrix discGram(GramMatrix(minimalBasis, args.partitions));
        outStream << "minBasisGramMatrix = "
            << MathematicaOutput(discGram) << std::endl;
        // DMatrix discGram_BasisStates = 
                // discPolys.transpose() * discGram * discPolys;
        // outStream << discGram << std::endl << "basis states:" 
                // << std::endl << discGram_BasisStates << std::endl;
    }

    timer.Start();
    DMatrix monoMassMatrix(MassMatrix(minimalBasis, args.partitions));
    DMatrix polyMassMatrix = discPolys.transpose()*monoMassMatrix*discPolys;

    if (outStream.rdbuf() != std::cout.rdbuf()) {
        outStream << "minBasisMassMatrix = "
            << MathematicaOutput(monoMassMatrix) << std::endl;
        outStream << "basisStateMassMatrix = "
                << MathematicaOutput(polyMassMatrix) << std::endl;
        std::cout << "Mass matrix computed in " << timer.TimeElapsedInWords()
            << "." << std::endl;
    } else {
        EigenSolver solver(polyMassMatrix.cast<builtin_class>());
        outStream << "Computed a mass matrix for the basis in " 
                << timer.TimeElapsedInWords() << "; its eigenvalues are:\n"
                << solver.eigenvalues() << std::endl;
    }

    timer.Start();
    DMatrix monoKineticMatrix(KineticMatrix(minimalBasis, args.partitions));
    DMatrix polyKineticMatrix = discPolys.transpose()*monoKineticMatrix*discPolys;

    if (outStream.rdbuf() != std::cout.rdbuf()) {
        outStream << "minBasisKineticMatrix = "
            << MathematicaOutput(monoKineticMatrix) << std::endl;
        outStream << "basisStateKineticMatrix = "
                << MathematicaOutput(polyKineticMatrix) << std::endl;
        std::cout << "Kinetic matrix computed in " << timer.TimeElapsedInWords()
            << "." << std::endl;
    } else {
        EigenSolver solver(polyKineticMatrix.cast<builtin_class>());
        outStream << "Computed a kinetic matrix for the basis in " 
                << timer.TimeElapsedInWords() << "; its eigenvalues are:\n"
                << solver.eigenvalues() << std::endl;
    }

    coeff_class m = 1;
    DMatrix hamiltonian = (m*m)*polyMassMatrix + polyKineticMatrix;
    if (outStream.rdbuf() != std::cout.rdbuf()) {
        outStream << "hamiltonian = "
                << MathematicaOutput(hamiltonian) << std::endl;
    } else {
        EigenSolver solver(hamiltonian.cast<builtin_class>());
        outStream << "Here are the Hamiltonian eigenvalues:\n" 
            << solver.eigenvalues() << std::endl;
    }


    return hamiltonian;
}

void GSLErrorHandler(const char* reason, const char* file, int line, int err) {
    std::cerr << "GSL Error in " << file << ":" << line << " --- "
        << gsl_strerror(err) << ", " << reason << std::endl;
}
