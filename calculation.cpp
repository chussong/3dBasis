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
    *args.outStream << "(*Hamiltonian test with delta=" << args.delta << "*)"
        << endl;
    
    *args.outStream << "(*EVEN STATES*)" << endl;
    Hamiltonian evenHam = FullHamiltonian(args, false);

    *args.outStream << "(*ODD STATES*)" << endl;
    Hamiltonian oddHam  = FullHamiltonian(args, true);

    *args.outStream << endl;

    return DMatrix();
}

// compute the hamiltonian for all states with delta up to args.delta; if 
// args.delta == 0, only compute one n-level (the DiagonalBlock at n=args.numP)
Hamiltonian FullHamiltonian(Arguments args, const bool odd) {
    int minN, maxN;
    if (args.delta != 0) {
        minN = 2;
        maxN = std::ceil(double(args.delta) / 1.5);
    } else {
        minN = args.numP;
        maxN = args.numP;
    }
    Hamiltonian output;
    output.maxN = maxN;
    const std::string parity = odd ? ", odd" : ", even";
    const bool mathematica = (args.options & OPT_MATHEMATICA) != 0;
    OStream& outStream = *args.outStream;

    std::vector<Basis<Mono>> minBases;
    std::vector<DMatrix> discPolys;
    for (int n = minN; n <= maxN; ++n) {
        // FIXME: remove adjustment so degree's consistently "L above dirichlet"
        if (args.delta != 0) {
            args.numP = n;
            args.degree = std::ceil(double(args.delta) - 0.5*n);
        } else {
            args.degree = args.degree + n;
        }

        // FIXME: directly generate only the monomials with the correct parity
        std::vector<Basis<Mono>> allEvenBases;
        std::vector<Basis<Mono>> allOddBases;
        for(int deg = n; deg <= args.degree; ++deg){
            splitBasis<Mono> degBasis(n, deg, args);
            allEvenBases.push_back(degBasis.EvenBasis());
            allOddBases.push_back(degBasis.OddBasis());
        }
        const std::vector<Basis<Mono>>& inputBases = 
                                            (odd ? allOddBases : allEvenBases);

        const std::string suffix = std::to_string(n) + parity;
        std::vector<Poly> orthogonalized = 
                        ComputeBasisStates_SameParity(inputBases, args, odd);
        minBases.push_back(MinimalBasis(orthogonalized));
        DMatrix polysOnMinBasis = PolysOnMinBasis(minBases[n-minN],
                                               orthogonalized, outStream);
        discPolys.push_back(DiscretizePolys(polysOnMinBasis, args.partitions));
        if (mathematica) {
            outStream << "minimalBasis[" << suffix << "] = "
                << MathematicaOutput(minBases[n-minN]) << endl;
            outStream << "(*Polynomials on this basis (as rows, not columns!):*)\n"
                << "polysOnMinBasis[" << suffix << "] = " 
                << MathematicaOutput(polysOnMinBasis.transpose()) << endl;
            outStream << "(*And discretized:*)\ndiscretePolys[" << suffix 
                << "] = " << MathematicaOutput(discPolys[n-minN].transpose()) 
                << endl;
        } else {
            outStream << "Minimal basis (" << n << "):" << minBases[n-minN] << endl;
        }

        output.diagonal.push_back(DiagonalBlock(minBases[n-minN], 
                                                discPolys[n-minN], 
                                                args, odd));
        if (n-2 >= minN) {
            output.nPlus2.push_back(NPlus2Block(minBases[n-2-minN], 
                                                discPolys[n-2-minN],
                                                minBases[n-minN],
                                                discPolys[n-minN], 
                                                args, odd));
        }
    }

    return output;
}

DMatrix DiagonalBlock(const Basis<Mono>& minimalBasis, 
                      const DMatrix& discPolys, 
                      const Arguments& args, const bool odd) {
    Timer timer;
    const bool interacting = (args.options & OPT_INTERACTING) != 0;
    std::string suffix = std::to_string(args.numP) + (odd ? ", odd" : ", even");

    timer.Start();
    DMatrix monoMassMatrix(MassMatrix(minimalBasis, args.partitions));
    DMatrix polyMassMatrix = discPolys.transpose()*monoMassMatrix*discPolys;
    OutputMatrix(monoMassMatrix, polyMassMatrix, "mass matrix", suffix, timer,
                 args);

    timer.Start();
    DMatrix monoKineticMatrix(KineticMatrix(minimalBasis, args.partitions));
    DMatrix polyKineticMatrix = discPolys.transpose()*monoKineticMatrix*discPolys;
    OutputMatrix(monoKineticMatrix, polyKineticMatrix, "kinetic matrix", suffix,
                 timer, args);

    DMatrix hamiltonian = args.msq*polyMassMatrix + polyKineticMatrix;
    if (interacting) {
        timer.Start();
        DMatrix monoNtoN(InteractionMatrix(minimalBasis, args.partitions));
        DMatrix polyNtoN = discPolys.transpose()*monoNtoN*discPolys;
        OutputMatrix(monoNtoN, polyNtoN, "NtoN matrix", suffix, timer, 
                     args);
        hamiltonian += args.lambda*polyNtoN;
    }

    /*
    if (mathematica) {
        outStream << "hamiltonian[" << suffix <<"] = "
                << MathematicaOutput(hamiltonian) << endl;
    } else {
        EigenSolver solver(hamiltonian.cast<builtin_class>());
        outStream << "Here are the Hamiltonian eigenvalues:\n" 
            << solver.eigenvalues() << endl;
    }
    */

    return hamiltonian;
}

// basisA is the minBasis of degree n, while basisB is the one for degree n+2
DMatrix NPlus2Block(const Basis<Mono>& basisA, const DMatrix& discPolysA,
                    const Basis<Mono>& basisB, const DMatrix& discPolysB,
                    const Arguments& args, const bool odd) {
    Timer timer;
    std::string suffix = std::to_string(args.numP) + (odd ? ", odd" : ", even");

    timer.Start();
    DMatrix monoNPlus2(NPlus2Matrix(basisA, basisB, args.partitions));
    DMatrix polyNPlus2 = discPolysA.transpose()*monoNPlus2*discPolysB;
    // FIXME?? if above line throws at runtime, swap A and B
    OutputMatrix(monoNPlus2, polyNPlus2, "NPlus2 matrix", suffix, timer, args);

    return args.lambda * polyNPlus2;
}

void OutputMatrix(const DMatrix& monoMatrix, const DMatrix& polyMatrix,
                  std::string name, const std::string& suffix, Timer& timer, 
                  const Arguments& args) {
    OStream& outStream = *args.outStream;
    OStream& console = *args.console;
    const bool mathematica = (args.options & OPT_MATHEMATICA) != 0;

    if (mathematica) {
        std::string mathematicaName = MathematicaName(name);
        name[0] = std::toupper(name[0]);
        outStream << "minBasis" << mathematicaName << "[" << suffix <<"] = "
            << MathematicaOutput(monoMatrix) << endl;
        outStream << "basisState" << mathematicaName << "[" << suffix <<"] = "
            << MathematicaOutput(polyMatrix) << endl;
        console << name << " computed in " << timer.TimeElapsedInWords()
            << "." << endl;
    } else {
        EigenSolver solver(polyMatrix.cast<builtin_class>());
        outStream << "Computed a " << name << " for the basis in " 
                << timer.TimeElapsedInWords() << "; its eigenvalues are:\n"
                << solver.eigenvalues() << endl;
    }
}

// capitalize each word, then delete all non-alphanumeric chars (inc. spaces)
std::string MathematicaName(std::string name) {
    name[0] = std::toupper(name[0]);
    bool capNext = false;
    for (std::size_t i = 0; i < name.size(); ++i) {
        if (capNext) {
            capNext = false;
            name[i] = std::toupper(name[i]);
        } else if (name[i] == ' ') {
            capNext = true;
        }
    }
    name.erase(std::remove_if(name.begin(), name.end(), 
                              [](char c){ return !std::isalnum(c); }),
               name.end());
    return name;
}

void GSLErrorHandler(const char* reason, const char* file, int line, int err) {
    std::cerr << "GSL Error in " << file << ":" << line << " --- "
        << gsl_strerror(err) << ", " << reason << std::endl;
}
