#include "calculation.hpp"

int Calculate(const Arguments& args) {
    gsl_set_error_handler(&GSLErrorHandler);

    if (args.options & OPT_TEST) {
        return Test::RunAllTests(args);
    }

    if (args.options & OPT_STATESONLY) {
        ComputeBasisStates(args);
        return EXIT_SUCCESS;
    }

    DMatrix hamiltonian = ComputeHamiltonian(args);
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

// output a matrix where each column is one of the basis vectors expressed in
// terms of the monomials on the minimal basis
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
    if (args.delta == 0.0) {
        *args.outStream << "(*Hamiltonian test at (n,l)=(" << args.numP 
            << "," << args.degree << "), ";
    } else {
        *args.outStream << "(*Hamiltonian test with delta=" << args.delta 
            << ", ";
    }
    *args.outStream << "kMax=" << args.partitions 
        << ". (m^2, \\lambda, \\Lambda) = (" << args.msq << ',' << args.lambda
        << ',' << args.cutoff << ")*)" << endl;

    Timer overallTimer;
    
    const bool full = (args.options & OPT_FULLOUTPUT) != 0;

    if (args.basisDir.empty()) {
        throw std::logic_error(__FILE__ ": ComputeHamiltonian called with empty"
                               " basisDir");
    }

    boost::filesystem::path basisPath = args.basisDir;

    *args.outStream << "(*EVEN STATES*)" << endl;
    Hamiltonian evenHam = FullHamiltonian(basisPath / "even", args, false);
    if (full) AnalyzeHamiltonian(evenHam, args, false);

    *args.outStream << "(*ODD STATES*)" << endl;
    Hamiltonian oddHam  = FullHamiltonian(basisPath / "odd", args, true);
    if (full) AnalyzeHamiltonian(oddHam, args, true);

    *args.console << "\nEntire computation took " 
        << overallTimer.TimeElapsedInWords() << "." << endl;

    return DMatrix();
}

namespace {
std::vector<Poly> NEqualsOneState() {
    std::vector<particle> basisParticle;
    basisParticle.push_back({1,0});

    Mono basisMono(basisParticle);
    builtin_class norm = InnerFock(basisMono, basisMono);
    basisMono /= std::sqrt(norm);

    std::vector<Poly> output;
    output.emplace_back(basisMono);

    return output;
}
} // anonymous namespace

// compute the hamiltonian for all states with delta up to args.delta; if 
// args.delta == 0, only compute one n-level (the DiagonalBlock at n=args.numP)
Hamiltonian FullHamiltonian(const boost::filesystem::path& basisDir, 
                            Arguments args, const bool odd) {
    Hamiltonian output;
    const std::string parity = odd ? ", odd" : ", even";
    const bool mathematica = (args.options & OPT_MATHEMATICA) != 0;
    OStream& outStream = *args.outStream;

    int highestN = -1;
    std::unordered_map<int,Basis<Mono>> minBases;
    std::unordered_map<int,SMatrix> discPolys;

    // do N=1 separately since it's trivial
    std::vector<Poly> nEqualsOne = NEqualsOneState();
    minBases.emplace(1, MinimalBasis(nEqualsOne));
    DMatrix polysOnMinBasis = PolysOnMinBasis(minBases.at(1), nEqualsOne,
                                              outStream);

    for (boost::filesystem::directory_entry& file 
            : boost::filesystem::directory_iterator(basisDir)) {
        int n = std::stoi(file.path().stem().string());
        highestN = std::max(highestN, n);
        std::vector<Poly> orthogonalized = Poly::ReadFromFile(file.path().string(), n);
        minBases.emplace(n, MinimalBasis(orthogonalized));
        DMatrix polysOnMinBasis = PolysOnMinBasis(minBases.at(n), orthogonalized,
                                                  outStream);

        if (n == 1) {
            SMatrix disc(polysOnMinBasis.rows(), polysOnMinBasis.cols());
            for (Eigen::Index row = 0; row < disc.rows(); ++row) {
                for (Eigen::Index col = 0; col < disc.cols(); ++col) {
                    disc.insert(row, col) = polysOnMinBasis(row, col);
                }
            }
            discPolys.emplace(n, std::move(disc));
        } else {
            discPolys.emplace(n, DiscretizePolys(polysOnMinBasis, args.partitions));
        }
        if (mathematica) {
            const std::string suffix = std::to_string(n) + parity;
            outStream << "minimalBasis[" << suffix << "] = "
                << MathematicaOutput(minBases.at(n)) << endl;
            outStream 
                << "(*Polynomials on this basis (as rows, not columns!):*)\n"
                << "polysOnMinBasis[" << suffix << "] = " 
                << MathematicaOutput(polysOnMinBasis.transpose()) << endl;
            outStream << "(*And discretized:*)\ndiscretePolys[" << suffix 
                << "] = " 
                << MathematicaOutput(DMatrix(discPolys.at(n).transpose())) 
                << endl;
        } else {
            outStream << "Minimal basis (" << n << "):" << minBases.at(n) 
                << endl;
        }

        args.numP = n;
        output.blocks.emplace(std::array<int,2>{n, n}, 
                              DiagonalBlock(minBases.at(n), discPolys.at(n), 
                                            args, odd));
        if ((args.options & OPT_INTERACTING) != 0) {
            if (discPolys.count(n-2) == 1) {
                args.numP = n;
                output.blocks.emplace(std::array<int,2>{n-2, n},
                                      NPlus2Block(minBases.at(n-2), discPolys.at(n-2),
                                                  minBases.at(n), discPolys.at(n),
                                                  args, odd));
            }
            if (discPolys.count(n+2) == 1) {
                args.numP = n+2;
                output.blocks.emplace(std::array<int,2>{n, n+2},
                                      NPlus2Block(minBases.at(n), discPolys.at(n),
                                                  minBases.at(n+2), discPolys.at(n+2),
                                                  args, odd));
            }
        }
    }

    // make sure hamiltonian has one block for each particle number up to the
    // highest, and assign startLocs to each block
    Eigen::Index runningCount = 0;
    for (int i = 1; i <= highestN; ++i) {
        if (output.blocks.count({i,i}) == 0) {
            throw std::runtime_error(__FILE__ ": Hamiltonian load incomplete: "
                                     "missing " + std::to_string(i) + ".txt");
        }

        output.startLocs.emplace(std::array<int,2>{i,i}, 
                                 std::array<Eigen::Index,2>{runningCount, runningCount});

        if (i-2 >= 1) {
            Eigen::Index trailingCount = runningCount - 
                                         output.blocks.at({i-2,i}).rows() -
                                         output.blocks.at({i-1,i}).rows();
            output.startLocs.emplace(std::array<int,2>{i-2,i}, 
                                     std::array<Eigen::Index,2>{trailingCount, runningCount});
        }

        runningCount += output.blocks.at({i,i}).rows();
    }

    return output;
}

DMatrix DiagonalBlock(const Basis<Mono>& minimalBasis, 
                      const SMatrix& discPolys, 
                      const Arguments& args, const bool odd) {
    *args.console << "DiagonalBlock(" << args.numP << ", " << args.degree << ")" 
        << endl;
    if (minimalBasis.size() == 0) return DMatrix(0, 0);

    Timer timer;
    const bool interacting = (args.options & OPT_INTERACTING) != 0;
    std::string suffix = std::to_string(args.numP) + (odd ? ", odd" : ", even");
    std::size_t kMax = (args.numP == 1 ? 1 : args.partitions);

    timer.Start();
    DMatrix monoMassMatrix(MassMatrix(minimalBasis, kMax));
    DMatrix polyMassMatrix = discPolys.transpose()*monoMassMatrix*discPolys;
    OutputMatrix(monoMassMatrix, polyMassMatrix, "mass matrix", suffix, timer,
                 args);

    timer.Start();
    DMatrix monoKineticMatrix(KineticMatrix(minimalBasis, kMax));
    DMatrix polyKineticMatrix = discPolys.transpose()*monoKineticMatrix*discPolys;
    OutputMatrix(monoKineticMatrix, polyKineticMatrix, "kinetic matrix", suffix,
                 timer, args);

    DMatrix hamiltonian = args.msq*polyMassMatrix 
                        + (args.cutoff*args.cutoff)*polyKineticMatrix;
    if (interacting) {
        timer.Start();
        DMatrix monoNtoN(InteractionMatrix(minimalBasis, kMax));
        DMatrix polyNtoN = discPolys.transpose()*monoNtoN*discPolys;
        OutputMatrix(monoNtoN, polyNtoN, "NtoN matrix", suffix, timer, 
                     args);
        hamiltonian += (args.lambda*args.cutoff)*polyNtoN;

        // if NtoN has nonpositive eigenvalues, it can't be right
        DEigenSolver solver(polyNtoN.cast<builtin_class>());
        if (solver.eigenvalues().rows() != 0 && solver.eigenvalues()(0) < 0
                && std::abs(solver.eigenvalues()(0)) > EPSILON) {
            *args.console << "Error: " << args.numP << "->" << args.numP 
                << " matrix not positive definite.\n";
        }
    }

    /*
    if (mathematica) {
        outStream << "hamiltonian[" << suffix <<"] = "
                << MathematicaOutput(hamiltonian) << endl;
    } else {
        DEigenSolver solver(hamiltonian.cast<builtin_class>());
        outStream << "Here are the Hamiltonian eigenvalues:\n" 
            << solver.eigenvalues() << endl;
    }
    */

    return hamiltonian;
}

// basisA is the minBasis of degree n, while basisB is the one for degree n+2
DMatrix NPlus2Block(const Basis<Mono>& basisA, const SMatrix& discPolysA,
                    const Basis<Mono>& basisB, const SMatrix& discPolysB,
                    const Arguments& args, const bool odd) {
    *args.console << "NPlus2Block(" << args.numP-2 << " -> " << args.numP << ")" 
        << endl;
    if (basisA.size() == 0 || basisB.size() == 0) return DMatrix(0, 0);

    Timer timer;
    std::string suffix = std::to_string(args.numP-2) 
                       + (odd ? ", odd" : ", even");

    timer.Start();
    DMatrix monoNPlus2(NPlus2Matrix(basisA, basisB, args.partitions));
    DMatrix polyNPlus2 = discPolysA.transpose()*monoNPlus2*discPolysB;
    OutputMatrix(monoNPlus2, polyNPlus2, "NPlus2 matrix", suffix, timer, args);

    return (args.lambda*args.cutoff) * polyNPlus2;
}

void AnalyzeHamiltonian(const Hamiltonian& hamiltonian, const Arguments& args,
                        const bool odd) {
    // Eigen::Index totalSize = 0;
    // for (const auto& block : hamiltonian.diagonal) totalSize += block.rows();
    // if (totalSize <= MAX_DENSE_SIZE) {
        // AnalyzeHamiltonian_Dense(hamiltonian, args);
    // } else {
        AnalyzeHamiltonian_Sparse(hamiltonian, args, odd);
    // }
}

// void AnalyzeHamiltonian_Dense(const Hamiltonian& hamiltonian, 
                               // const Arguments& args) {
    // Eigen::Index totalSize = 0;
    // for (const auto& block : hamiltonian.diagonal) totalSize += block.rows();
    // DMatrix matrixForm(totalSize, totalSize);
// 
    // Eigen::Index offset = 0;
    // Eigen::Index trailingOffset = 0;
    // for (std::size_t n = 1; n < hamiltonian.diagonal.size()+1; ++n) {
        // const auto& block = hamiltonian.diagonal[n-1];
        // matrixForm.block(offset, offset, block.rows(), block.cols()) = block;
// 
        // if (n >= 3) {
            // const auto& nPlus2Block = hamiltonian.nPlus2[n-4];
            // matrixForm.block(trailingOffset, offset, nPlus2Block.rows(),
                             // nPlus2Block.cols()) = nPlus2Block;
            // matrixForm.block(offset, trailingOffset, nPlus2Block.cols(),
                             // nPlus2Block.rows()) = nPlus2Block.transpose();
            // trailingOffset += nPlus2Block.rows();
        // }
        // offset += block.rows();
    // }
// 
    // DEigenSolver solver(matrixForm.cast<builtin_class>());
    // *args.console << "Hamiltonian eigenvalues:\n" 
        // << solver.eigenvalues() << endl;
// }

void AnalyzeHamiltonian_Sparse(const Hamiltonian& hamiltonian, 
                               const Arguments& args, const bool odd) {
    Eigen::Index maxIndex = -1;
    std::vector<Triplet> triplets;
    for (const auto& keyValuePair : hamiltonian.blocks) {
        // get which block this is and the location of its top left corner
        const DMatrix& block = keyValuePair.second;
        const std::array<Eigen::Index,2>& startLoc = 
            hamiltonian.startLocs.at(keyValuePair.first);

        // make sure to get the maximum size of the matrix
        Eigen::Index extentOfBlock = std::max(startLoc[0] + block.rows(),
                                              startLoc[1] + block.cols());
        maxIndex = std::max(maxIndex, extentOfBlock);

        for (Eigen::Index i = 0; i < block.rows(); ++i) {
            for (Eigen::Index j = 0; j < block.cols(); ++j) {
                if (block(i,j) == 0) continue;

                triplets.emplace_back(startLoc[0]+i, startLoc[1]+j, block(i,j));
                // std::cout << "Triplet: " << triplets.back() << '\n';
                if (startLoc[0] != startLoc[1]) {
                    triplets.emplace_back(startLoc[1]+j, startLoc[0]+i, 
                                          block(i,j));
                }
                // std::cout << "Triplet: " << triplets.back() << '\n';
            }
        }
    }

    SMatrix matrixForm(maxIndex, maxIndex);
    // std::cout << "Matrix size: " << maxIndex << "x" << maxIndex << std::endl;
    // for (const auto& trip : triplets) std::cout << trip << std::endl;
    matrixForm.setFromTriplets(triplets.begin(), triplets.end());

    OStream& outStream = *args.outStream;
    if ((args.options & OPT_MATHEMATICA) != 0) {
        outStream << "hamiltonian[" << (odd ? "odd" : "even") << "] = "
            << MathematicaOutput(matrixForm) << '\n';
    } else {
        if (matrixForm.rows() <= 20) {
            outStream << (odd ? "Odd" : "Even") << " Hamiltonian:\n" 
                << DMatrix(matrixForm) << '\n';
        }
        outStream << (odd ? "Odd" : "Even") << " Hamiltonian eigenvalues:\n";
        // TODO: use sparse Lanczos instead of converting to dense?
        DMatrix denseForm(matrixForm);
        DEigenSolver solver(denseForm.cast<builtin_class>());
        outStream << solver.eigenvalues() << '\n';
    }
}

void OutputMatrix(const DMatrix& monoMatrix, const DMatrix& polyMatrix,
                  std::string name, const std::string& suffix, Timer& timer, 
                  const Arguments& args) {
    OStream& outStream = *args.outStream;
    OStream& console = *args.console;
    const bool mathematica = (args.options & OPT_MATHEMATICA) != 0;
    const bool full = (args.options & OPT_FULLOUTPUT) != 0;

    if (mathematica) {
        std::string mathematicaName = MathematicaName(name);
        name[0] = std::toupper(name[0]);
        if (full) {
            outStream << "minBasis" << mathematicaName << "[" << suffix <<"] = "
                << MathematicaOutput(monoMatrix) << '\n';
        }
        outStream << "basisState" << mathematicaName << "[" << suffix << "] = "
            << MathematicaOutput(polyMatrix) << '\n';
        console << name << " computed in " << timer.TimeElapsedInWords()
            << ".\n";
        // if (polyMatrix.rows() == polyMatrix.cols()) {
            // DEigenSolver solver(polyMatrix.cast<builtin_class>());
            // console << name << " computed in " << timer.TimeElapsedInWords()
                // << "; its eigenvalues are:\n" << solver.eigenvalues() << '\n';
        // }
    } else if (polyMatrix.rows() <= 10 && polyMatrix.cols() <= 10) {
        outStream << "Computed a " << name << " for the basis in " 
            << timer.TimeElapsedInWords() << "; ";
        if (full) outStream << "mono:\n" << monoMatrix << '\n';
        outStream << "poly:\n" << polyMatrix << '\n';
    } else if (polyMatrix.rows() == polyMatrix.cols()) {
        DEigenSolver solver(polyMatrix.cast<builtin_class>());
        outStream << "Computed a " << name << " for the basis in " 
            << timer.TimeElapsedInWords() << "; its eigenvalues are:\n"
            << solver.eigenvalues() << '\n';
    } else {
        outStream << "Computed a " << name << " for the basis in " 
            << timer.TimeElapsedInWords() << ", but it's not square and is "
            << "too large to show." << '\n';
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
