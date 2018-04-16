#include "main.hpp"

int main(int argc, char* argv[]) {
    Arguments args = ParseArguments(argc, argv);
    if (args.options & OPT_VERSION){
        std::cout << "This is 3dBasis version " << VERSION << ", released "
            << RELEASE_DATE << ". The latest updates can always be found at "
            << "https://github.com/chussong/3dBasis." << std::endl;
        return EXIT_SUCCESS;
    }

    if (args.degree == 0 || args.numP == 0) {
        // std::cerr << "Error: you must enter a number of particles and a degree."
            // << std::endl;
        // return EXIT_FAILURE;
        return GUI::StartGUI(argc, argv, args);
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
