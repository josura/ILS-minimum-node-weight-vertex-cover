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

        uint removesDepth;

        bool recursiveDepth;

    public:
        //TESTING
        uint globalRemove = 0,globalSwap = 0,globalHybrid = 0; 
        //TESTING

        LocalSearch(WeightedVertexGraph* _graph,double samplingFactor=7,int depthOfRemoves = 50, bool recursiveDepth = true, uint numberOfIterations=500);

        ~LocalSearch();     

        NodeBitArray startResolve();
        NodeBitArray startResolveStrange();
        NodeBitArray startResolveOptimized();
        NodeBitArray startResolveOptimized(double &finalCost,NodeBitArray startSolution=nullptr);
        NodeBitArray startResolveWithLimit(double &finalCost,NodeBitArray startSolution=nullptr);
        NodeBitArray startResolveWithLimitAndSampling(double &finalCost,NodeBitArray startSolution=nullptr);
        NodeBitArray startResolveWithLimitAndSamplingOptimized(double &finalCost,NodeBitArray startSolution=nullptr);
        NodeBitArray startResolveFinal(double &finalCost,NodeBitArray startSolution=nullptr);

        void setSolution(NodeBitArray solution);
        
        NodeBitArray getSolution()const;
        NodeSet getSolutionSet()const;

        double getSolutionWeight();

        uint getCurrentNumberOfObjectiveFunctionEval()const;
};