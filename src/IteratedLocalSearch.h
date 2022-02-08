#pragma once 


#include "WeightedVertexGraph.h"
#include "LocalSearch.h"
#include "utilities.h"
#include "heuristics.h"
#include <sys/types.h>

class IteratedLocalSearch{
    private:
        NodeBitArray solution;

        //TODO see if memorizing the best solution found is the best thing or not
        NodeBitArray currentBestSolution;

        double solutionCost;

        //parameters
        double power;

        double rarity;

        double probabilityBadSolution;

        uint maxNumberOfIterations;

        uint maxNumberOfIterationsWithouthImprovement;

        WeightedVertexGraph* graph;

        LocalSearch* localSearchMethod;

        NodeBitArray bestSolution;



    public:
        IteratedLocalSearch(WeightedVertexGraph* _graph,uint numberOfIterations=500, uint numberofIterationsWithouthImprovement=5, 
            double startingPower=10, 
            double startingRarity = 50, 
            double startingProbabilityBadSolution=5
        );

        ~IteratedLocalSearch();     

        NodeBitArray simplePerturbationSolution(NodeBitArray solution);

        double adjustPerturbationPowerAndRarity();


        double adjustProbabilityOfBadSolutions();

        NodeBitArray startResolve();
        
        NodeBitArray getSolution()const;
        NodeSet* getSolutionSet()const;

        double getSolutionWeight()const;
};