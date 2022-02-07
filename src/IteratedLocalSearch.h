#pragma once 


#include "WeightedVertexGraph.h"
#include "LocalSearch.h"
#include "utilities.h"
#include "heuristics.h"

class IteratedLocalSearch{
    private:
        NodeBitArray solution;

        NodeSet* solutionSet;

        uint numberOfIterations;

        WeightedVertexGraph* graph;

    public:
        IteratedLocalSearch(WeightedVertexGraph* _graph,uint numberOfIterations=500);

        ~IteratedLocalSearch();     

        NodeBitArray startResolve();
        
        NodeBitArray getSolution()const;
        NodeSet* getSolutionSet()const;

        double getSolutionWeight()const;
};