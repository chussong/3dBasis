#include "calculation.hpp"

int Calculate(const Arguments& args) {
    OStream& console = *args.console;
    gsl_set_error_handler(&GSLErrorHandler);

    if (args.options & OPT_TEST) {
        return Test::RunAllTests(args);
    }

    // initialize all multinomials which might come up
    //
    // this is obviously something of a blunt instrument and could easily be
    // made more efficient

    for (int n = 1; n <= args.numP; ++n) {
        console << "Initialize(" << n << ", " << 2*args.degree << ")" 
            << endl;
        Multinomial::Initialize(n, 2*args.degree);
    }

    if (args.options & OPT_STATESONLY) {
        ComputeBasisStates(args);
        // if (args.outStream->rdbuf() != std::cout.rdbuf()) {
            // delete args.outStream;
        // }
        return EXIT_SUCCESS;
    }

    DMatrix hamiltonian = ComputeHamiltonian(args);
    // if (args.outStream->rdbuf() != std::cout.rdbuf()) {
        // delete args.outStream;
    // }
    return EXIT_SUCCESS;
}

// return basis polynomials. They are NOT normalized w.r.t. partitions
std::vector<Poly> ComputeBasisStates(const Arguments& args) {
    int numP = args.numP;
    int degree = args.degree + args.numP; // add required Dirichlet derivatives
    // int options = args.options;

    *args.outStream << "(*Orthogonal basis states with N=" << numP << ", L="
        << degree << " (including Dirichlet derivatives).*)" << endl;
    
    std::vector<Basis<Mono>> allEvenBases;
    std::vector<Basis<Mono>> allOddBases;
    for (int deg = numP; deg <= degree; ++deg) {
        splitBasis<Mono> degBasis(numP, deg, args);
        allEvenBases.push_back(degBasis.EvenBasis());
        allOddBases.push_back(degBasis.OddBasis());
    }

    *args.outStream << "(*EVEN STATE ORTHOGONALIZATION*)" << endl;
    std::vector<Poly> basisEven = ComputeBasisStates_SameParity(allEvenBases, 
                                                                args, false);

    *args.outStream << "(*ODD STATE ORTHOGONALIZATION*)" << endl;
    std::vector<Poly> basisOdd = ComputeBasisStates_SameParity(allOddBases, 
                                                               args, true);

    *args.outStream << endl;

    basisEven.insert(basisEven.end(), basisOdd.begin(), basisOdd.end());
    return basisEven;
}

// return basis polynomials. They are NOT normalized w.r.t. partitions
std::vector<Poly> ComputeBasisStates_SameParity(
        const std::vector<Basis<Mono>>& inputBases, const Arguments& args,
        const bool odd) {
    OStream& console = *args.console;
    std::vector<Poly> orthogonalized = Orthogonalize(inputBases, console, odd);

    // Basis<Mono> minimalBasis(MinimalBasis(orthogonalized));
    // if (outStream.rdbuf() == std::cout.rdbuf()) {
        // outStream << "Minimal basis: " << minimalBasis << endl;
    // } else {
        // outStream << "minimalBasis[" << suffix <<"] = " << MathematicaOutput(minimalBasis) 
            // << endl;
    // }

    return orthogonalized;
}

DMatrix PolysOnMinBasis(const Basis<Mono>& minimalBasis,
        const std::vector<Poly> orthogonalized, OStream&) {
    DMatrix polysOnMinBasis(minimalBasis.size(), orthogonalized.size());
    for (std::size_t i = 0; i < orthogonalized.size(); ++i) {
        polysOnMinBasis.col(i) = minimalBasis.DenseExpressPoly(
                        orthogonalized[i] );
    }

    // outStream << "(*Polynomials on this basis (as rows, not columns!):*)\n"
        // "polysOnMinBasis[" << suffix <<"] = " << MathematicaOutput(polysOnMinBasis.transpose()) 
        // << endl;

    return polysOnMinBasis;
}
    

DMatrix ComputeHamiltonian(const Arguments& args) {
    int numP = args.numP;
    int degree = args.degree + args.numP; // add required Dirichlet derivatives
    // int options = args.options;

    *args.outStream << "(*Matrix element test with N=" << numP << ", L="
        << degree << " (including Dirichlet derivatives).*)" << endl;
    
    std::vector<Basis<Mono>> allEvenBases;
    std::vector<Basis<Mono>> allOddBases;
    for(int deg = numP; deg <= degree; ++deg){
        splitBasis<Mono> degBasis(numP, deg, args);
        allEvenBases.push_back(degBasis.EvenBasis());
        allOddBases.push_back(degBasis.OddBasis());
    }

    *args.outStream << "(*EVEN STATE ORTHOGONALIZATION*)" << endl;
    DMatrix evenHam = ComputeHamiltonian_SameParity(allEvenBases, args, false);

    *args.outStream << "(*ODD STATE ORTHOGONALIZATION*)" << endl;
    DMatrix oddHam  = ComputeHamiltonian_SameParity(allOddBases, args, true);

    *args.outStream << endl;

    return DMatrix();
}

// the portion of the hamiltonian computation which assumes equal pt parity
DMatrix ComputeHamiltonian_SameParity(const std::vector<Basis<Mono>>& inputBases,
                                      const Arguments& args, const bool odd) {
    Timer timer;
    OStream& outStream = *args.outStream;
    OStream& console = *args.console;
    bool mathematica = (args.options & OPT_MATHEMATICA) != 0;
    bool interacting = (args.options & OPT_INTERACTING) != 0;

    std::vector<Poly> orthogonalized = ComputeBasisStates_SameParity(inputBases,
            args, odd);
    const std::string suffix = odd ? "odd" : "even";

    Basis<Mono> minimalBasis(MinimalBasis(orthogonalized));
    if (mathematica) {
        outStream << "minimalBasis[" << suffix << "] = " << 
            MathematicaOutput(minimalBasis) << endl;
    } else {
        outStream << "Minimal basis: " << minimalBasis << endl;
    }

    DMatrix polysOnMinBasis = PolysOnMinBasis(minimalBasis, orthogonalized,
            outStream);

    // outStream << "polysOnMinBasis:\n" << polysOnMinBasis << endl;
    DMatrix discPolys = DiscretizePolys(polysOnMinBasis, args.partitions);
    if (mathematica) {
        outStream << "(*Polynomials on this basis (as rows, not columns!):*)\n"
                << "polysOnMinBasis[" << suffix << "] = " 
        << MathematicaOutput(polysOnMinBasis.transpose()) << endl;
        outStream << "(*And discretized:*)\ndiscretePolys[" << suffix << "] = "
            << MathematicaOutput(discPolys.transpose()) << endl;
    }

    if (mathematica) {
        DMatrix discGram(GramMatrix(minimalBasis, args.partitions));
        outStream << "minBasisGramMatrix[" << suffix << "] = "
            << MathematicaOutput(discGram) << endl;
        // DMatrix discGram_BasisStates = 
                // discPolys.transpose() * discGram * discPolys;
        // outStream << discGram << endl << "basis states:" 
                // << endl << discGram_BasisStates << endl;
    }

    timer.Start();
    DMatrix monoMassMatrix(MassMatrix(minimalBasis, args.partitions));
    DMatrix polyMassMatrix = discPolys.transpose()*monoMassMatrix*discPolys;

    if (mathematica) {
        outStream << "minBasisMassMatrix[" << suffix <<"] = "
            << MathematicaOutput(monoMassMatrix) << endl;
        outStream << "basisStateMassMatrix[" << suffix <<"] = "
                << MathematicaOutput(polyMassMatrix) << endl;
        console << "Mass matrix computed in " << timer.TimeElapsedInWords()
            << "." << endl;
    } else {
        EigenSolver solver(polyMassMatrix.cast<builtin_class>());
        outStream << "Computed a mass matrix for the basis in " 
                << timer.TimeElapsedInWords() << "; its eigenvalues are:\n"
                << solver.eigenvalues() << endl;
    }

    timer.Start();
    DMatrix monoKineticMatrix(KineticMatrix(minimalBasis, args.partitions));
    DMatrix polyKineticMatrix = discPolys.transpose()*monoKineticMatrix*discPolys;

    if (mathematica) {
        outStream << "minBasisKineticMatrix[" << suffix <<"] = "
            << MathematicaOutput(monoKineticMatrix) << endl;
        outStream << "basisStateKineticMatrix[" << suffix <<"] = "
                << MathematicaOutput(polyKineticMatrix) << endl;
        console << "Kinetic matrix computed in " << timer.TimeElapsedInWords()
            << "." << endl;
    } else {
        EigenSolver solver(polyKineticMatrix.cast<builtin_class>());
        outStream << "Computed a kinetic matrix for the basis in " 
                << timer.TimeElapsedInWords() << "; its eigenvalues are:\n"
                << solver.eigenvalues() << endl;
    }

    DMatrix hamiltonian = args.msq*polyMassMatrix + polyKineticMatrix;
    if (interacting) {
        timer.Start();
        DMatrix monoNtoN(InteractionMatrix(minimalBasis, args.partitions));
        DMatrix polyNtoN = discPolys.transpose()*monoNtoN*discPolys;
        if (mathematica) {
            outStream << "minBasisNtoNMatrix[" << suffix <<"] = "
                << MathematicaOutput(monoNtoN) << endl;
            outStream << "basisStateNtoNMatrix[" << suffix <<"] = "
                    << MathematicaOutput(polyNtoN) << endl;
            console << "Interaction matrix computed in " 
                << timer.TimeElapsedInWords() << "." << endl;
        } else {
            outStream << "Computed an n-to-n interaction matrix for the basis "
                << "in " << timer.TimeElapsedInWords() << "." << endl;
        }

        timer.Start();
        DMatrix monoNPlus2(NPlus2Matrix(minimalBasis, args.partitions));
        DMatrix polyNPlus2 = discPolys.transpose()*monoNPlus2*discPolys;
        if (mathematica) {
            outStream << "minBasisNPlus2Matrix[" << suffix <<"] = "
                << MathematicaOutput(monoNPlus2) << endl;
            outStream << "basisStateNPlus2Matrix[" << suffix <<"] = "
                    << MathematicaOutput(polyNPlus2) << endl;
            console << "N+2 interaction matrix computed in " 
                << timer.TimeElapsedInWords() << "." << endl;
        } else {
            outStream << "Computed an n-to-n+2 interaction matrix for the "
                << "basis in " << timer.TimeElapsedInWords() << "." << endl;
        }

        hamiltonian += args.lambda*(polyNtoN + polyNPlus2);
    }

    if (mathematica) {
        outStream << "hamiltonian[" << suffix <<"] = "
                << MathematicaOutput(hamiltonian) << endl;
    } else {
        EigenSolver solver(hamiltonian.cast<builtin_class>());
        outStream << "Here are the Hamiltonian eigenvalues:\n" 
            << solver.eigenvalues() << endl;
    }

    return hamiltonian;
}

void GSLErrorHandler(const char* reason, const char* file, int line, int err) {
    std::cerr << "GSL Error in " << file << ":" << line << " --- "
        << gsl_strerror(err) << ", " << reason << std::endl;
}
