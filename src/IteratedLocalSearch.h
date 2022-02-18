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

        double propagationProportion;

        uint maxNumberOfIterations;

        uint maxNumberOfIterationsWithouthImprovement;

        WeightedVertexGraph* graph;

        LocalSearch* localSearchMethod;

        NodeBitArray bestSolution;

        //things required for the relation
        int objectiveEvalForOptimum=0;



    public:
        IteratedLocalSearch(WeightedVertexGraph* _graph,double startingPropagationProportion=5,uint numberOfIterations=100000, uint numberofIterationsWithouthImprovement=100000, 
            double startingPower=5, 
            double startingRarity = 50, 
            double startingProbabilityBadSolution=0.1
        );

        ~IteratedLocalSearch();     

        NodeBitArray simplePerturbationSolution(NodeBitArray solution);
        NodeBitArray randomPerturbationSolution(NodeBitArray solution);
        NodeBitArray randomDistributedPerturbationSolution(NodeBitArray solution);
        NodeBitArray randomDistributedPerturbationAndRemoves(NodeBitArray solution);

        void adjustPerturbationPowerAndRarity(int iterationWithoutImprovement, int iterationWithImprovement);


        void adjustProbabilityOfBadSolutions(int iterationWithoutImprovement, int iterationWithImprovement,int badSolutionsAccepted);

        NodeBitArray startResolve();
        
        NodeBitArray getSolution()const;
        NodeSet* getSolutionSet()const;

        double getSolutionWeight()const;

        double getPower()const;

        double getRarity()const;

        double getProbabilityBadSolution()const;

        int getObjectiveEvalForOptimum()const;
};