#include "main.hpp"

int main(int argc, char* argv[]) {
    Arguments args = ParseArguments(argc, argv);
    if (args.options & OPT_VERSION){
        *args.console << "This is 3dBasis version " << VERSION << ", released "
            << RELEASE_DATE << ". The latest updates can always be found at "
            << "https://github.com/chussong/3dBasis." << endl;
        return EXIT_SUCCESS;
    }

    if (args.degree == 0 || args.numP == 0) {
#ifdef NO_GUI
        std::cerr << "Error: you must enter a number of particles and a degree."
            << std::endl;
        return EXIT_FAILURE;
#else
        return GUI::StartGUI(argc, argv, args);
#endif
    }

    /* FIXME: move all of this to test.cpp
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
    */

    return Calculate(args);
}

Arguments ParseArguments(int argc, char* argv[]) {
    std::vector<std::string> options;
    std::string arg;
    Arguments ret;
    int j = 0;
    for (int i = 1; i < argc; ++i) {
        arg = argv[i];
        if (arg.size() > 0) {
            if (arg[0] == '-') {
                if (arg.size() > 1 && arg[1] == 'o') {
                    // open next argument as outstream, appending to it
#ifdef NO_GUI
                    ret.outStream = new std::ofstream(argv[i+1], 
                                    std::ios_base::out | std::ios_base::app);
#else
                    QFile* file = new QFile(argv[i+1]);
                    if (!file->open(QFile::WriteOnly | QFile::Append)) {
                        throw std::runtime_error("couldn't open file");
                    }
                    ret.outStream = new QTextStream(file);
#endif
                    ret.options |= OPT_MATHEMATICA;
                    ++i; // next argument is the filename so don't process it
                } else if (arg.size() > 1 && arg[1] == 'O') {
                    // open next argument as outstream, replacing it
#ifdef NO_GUI
                    ret.outStream = new std::ofstream(argv[i+1], 
                                    std::ios_base::out | std::ios_base::app);
#else
                    QFile* file = new QFile(argv[i+1]);
                    if (!file->open(QFile::WriteOnly | QFile::Truncate)) {
                        throw std::runtime_error("couldn't open file");
                    }
                    ret.outStream = new QTextStream(file);
#endif
                    ret.options |= OPT_MATHEMATICA;
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
                            ret.partitions = ReadArg<int>(arg);
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
#ifdef NO_GUI
    ret.console = new std::ostream(std::cout.rdbuf());
#else
    ret.console = new QTextStream(stdout);
#endif
    if (ret.outStream == nullptr) ret.outStream = ret.console;
    if(j < 2) ret.numP = 0; // invalidate the input since it was insufficient
    ret.options |= ParseOptions(options);
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
            ret |= OPT_INTERACTING;
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
