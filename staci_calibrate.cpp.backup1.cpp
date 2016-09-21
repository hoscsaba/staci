//
// Created by Hős Csaba on 2016. 09. 19..
//

#include <stdio.h>
#include <ga/ga.h>
#include <vector>
#include "Staci.h"

using namespace std;
#define cout STD_COUT

float Objective(GAGenome &);

Staci *feladat;
vector<string> pipe_names;

int main(int argc, char **argv) {

    pipe_names.push_back("PIPE18");
    pipe_names.push_back("PIPE4");

    feladat = new Staci(argc, argv);
    feladat->build_system();
    feladat->ini();
    bool success = feladat->solve_system();

    vector<double> origD;
    for (unsigned int i = 0; i < pipe_names.size(); i++)
        origD.push_back(feladat->get_dprop(pipe_names.at(i), "diameter"));

    // See if we've been given a seed to use (for testing purposes).  When you
    // specify a random seed, the evolution will be exactly the same each time
    // you use that seed number.

    unsigned int seed = 0;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i++], "seed") == 0) {
            seed = atoi(argv[i]);
        }
    }

    // Declare variables for the GA parameters and set them to some default values.

    int popsize = 30;
    int ngen = 100;
    float pmut = 0.01;
    float pcross = 0.6;

    // Create a phenotype for two variables.  The number of bits you can use to
    // represent any number is limited by the type of computer you are using.  In
    // this case, we use 16 bits to represent a floating point number whose value
    // can range from -5 to 5, inclusive.  The bounds on x1 and x2 can be applied
    // here and/or in the objective function.
    GABin2DecPhenotype map;
    double LOWER_BOUND,UPPER_BOUND;
    for (unsigned int i = 0; i < pipe_names.size(); i++) {
        LOWER_BOUND = 0.5 * feladat->get_dprop(pipe_names.at(i), "diameter");
        UPPER_BOUND = 2.0 * feladat->get_dprop(pipe_names.at(i), "diameter");
        map.add(16, LOWER_BOUND, UPPER_BOUND);
    }

    // Create the template genome using the phenotype map we just made.

    GABin2DecGenome genome(map, Objective);

    // Now create the GA using the genome and run it.  We'll use sigma truncation
    // scaling so that we can handle negative objective scores.

    GASimpleGA ga(genome);
    GASigmaTruncationScaling scaling;
    ga.populationSize(popsize);
    ga.nGenerations(ngen);
    ga.pMutation(pmut);
    ga.pCrossover(pcross);
    ga.scaling(scaling);
    ga.scoreFilename("bog.dat");
    ga.scoreFrequency(10);
    ga.flushFrequency(50);
    ga.evolve(seed);

    // Dump the results of the GA to the screen.

    genome = ga.statistics().bestIndividual();
    cout << "Solution found by the GA:";
    for (unsigned int i = 0; i < pipe_names.size(); i++) {
        printf("\n %10s: %5.3f (=? %5.3f, origD: %5.3f)",pipe_names.at(i).c_str(),feladat->get_dprop(pipe_names.at(i), "diameter"),genome.phenotype(i),origD.at(i));
    }
    cout << "best of generation data are in '" << ga.scoreFilename() << "'\n";

    return 0;
}

float
Objective(GAGenome &x) {

    GABin2DecGenome &genome = (GABin2DecGenome &) x;
    for (unsigned int i = 0; i < genome.nPhenotypes(); i++)
        feladat->set_dprop(pipe_names[i], "diameter", genome.phenotype(i));

    bool success = feladat->solve_system();

    float obj = 0.;
    if (success) {
        float mp = 0;
        for (unsigned int i = 0; i < genome.nPhenotypes(); i++) {
            mp = feladat->get_dprop(pipe_names.at(i).c_str(), "mass_flow_rate");
            obj += fabs(mp);
            printf("\n mass flow rate of %9s (D=%5.3fm) is %g, overall: %g",
                   pipe_names.at(i).c_str(), feladat->get_dprop(pipe_names[i], "diameter"), mp, obj);
        }
    } else{
        obj=1.e6;
        printf("\n\n THE SOLVER DID NOT CONVERGE!\n\n");
    }

    return obj;
}


