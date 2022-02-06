#pragma once 


#include "WeightedVertexGraph.h"
#include "utilities.h"
#include "heuristics.h"

class LocalSearch{
    private:
        NodeBitArray solution;

        uint numberOfIterations;

        WeightedVertexGraph* graph;

    public:
        LocalSearch(WeightedVertexGraph* _graph,uint numberOfIterations=500);

        ~LocalSearch();     

        NodeBitArray startResolveOptimized();
        
        NodeBitArray getSolution()const;

        double getSolutionWeight()const;
};