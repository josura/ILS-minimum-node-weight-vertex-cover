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

        uint currentObjectiveFunEvalution;

        uint maximumObjectiveFunEvalution = 20000;

        double samplingFactor;

    public:
        LocalSearch(WeightedVertexGraph* _graph,uint numberOfIterations=500,double samplingFactor=100);

        ~LocalSearch();     

        NodeBitArray startResolve();
        NodeBitArray startResolveStrange();
        NodeBitArray startResolveOptimized();
        NodeBitArray startResolveOptimized(double &finalCost,NodeBitArray startSolution=nullptr);
        NodeBitArray startResolveWithLimit(double &finalCost,NodeBitArray startSolution=nullptr);
        NodeBitArray startResolveWithLimitAndSampling(double &finalCost,NodeBitArray startSolution=nullptr);

        void setSolution(NodeBitArray solution);
        
        NodeBitArray getSolution()const;
        NodeSet getSolutionSet()const;

        double getSolutionWeight();

        uint getCurrentNumberOfObjectiveFunctionEval()const;
};