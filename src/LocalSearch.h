#pragma once 


#include "WeightedVertexGraph.h"
#include "utilities.h"
#include "heuristics.h"

class LocalSearch{
    private:
        NodeBitArray solution;

        NodeSet solutionSet;

        uint numberOfIterations;

        WeightedVertexGraph* graph;

    public:
        LocalSearch(WeightedVertexGraph* _graph,uint numberOfIterations=500);

        ~LocalSearch();     

        NodeBitArray startResolve();
        NodeBitArray startResolveStrange();
        NodeBitArray startResolveOptimized();
        NodeBitArray startResolveOptimized(double &finalCost,NodeBitArray startSolution=nullptr);

        void setSolution(NodeBitArray solution);
        
        NodeBitArray getSolution()const;
        NodeSet getSolutionSet()const;

        double getSolutionWeight()const;
};