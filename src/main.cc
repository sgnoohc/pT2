#include <exception>
#include <iostream>

#include "config.h"
#include "gator.h"
#include "plot.h"
#include "process.h"
#include "scan.h"

int main(int argc, char **argv)
{
    Config cfg;
    if (!parseArgs(argc, argv, cfg)) return 1;

    print_gator();
    printConfig(cfg);

    try
    {
        switch (cfg.mode)
        {
            case Mode::Process:
                return runProcess(cfg);
            case Mode::Scan:
                return runScan(cfg);
            case Mode::Plot:
                return runPlot(cfg);
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 1;
}
