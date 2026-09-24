#ifndef GATOR_H
#define GATOR_H

#include "RtypesCore.h"
#include <iostream>
#include <iomanip>
#include <cmath>

void print_creature();
void print_gator();

// Progress Bar

#define RED     "\033[1;31m"
#define YELLOW  "\033[1;33m"
#define GREEN   "\033[1;32m"
#define RESET   "\033[0m"

inline void printProgressBar(Long64_t current, Long64_t total)
{
    const int barWidth = 50;
    double progress = (double)current / (double)total;
    int pos = (int)(barWidth * progress);
    // Choose color
    const char* color;
    if (progress < 0.33)
        color = RED;
    else if (progress < 0.66)
        color = YELLOW;
    else
        color = GREEN;
    std::cout << "\r[";
    // Filled portion
    std::cout << color;
    for (int i = 0; i < pos; ++i)
        std::cout << "=";
    std::cout << RESET;
    // Empty portion
    for (int i = pos; i < barWidth; ++i)
        std::cout << " ";
    std::cout << "] ";
    std::cout << std::fixed << std::setprecision(1)
              << (progress * 100.0) << "%";
    std::cout.flush();
    // Finish line when done
    if (current == total - 1)
        std::cout << std::endl;
}

#endif
