#pragma once

#include <ostream>
#include <string>
#include <sys/types.h>
#include <unordered_set>
#include <iostream>
#include <vector>
#include <utility>

class WeightedVertexGraph{
    private:
        uint numberOfNodes;
        uint numberOfEdges=0;
        double* nodeWeights;
        std::unordered_set<uint>* adjList;
        std::vector<uint>* adjVector;
        std::vector<std::pair<uint, uint> > edgesVector;
        std::pair<uint, uint>* edgesArray;

    public:
        WeightedVertexGraph();

        WeightedVertexGraph(uint numNodes, double* nodeWeights);

        ~WeightedVertexGraph();

        WeightedVertexGraph* addEdge(uint node1, uint node2);

        std::pair<uint, uint>* makeEdgesArray();

        // accessory functions

        uint getNumNodes()const ;
        uint getNumEdges()const ;
        uint degreeOfNode(uint node)const;

        double* getNodeWeights()const;

        std::string getNodeWeightsStr()const;

        std::unordered_set<uint>* getAdjList(uint node)const;

        std::string getAdjListStr(uint node)const;

        bool adjNodes(uint node1, uint node2);

        std::vector<std::pair<uint, uint>> getEdgesVector()const;

        std::pair<uint, uint>* getEdgesArray()const;

        double getNodeWeight(uint node)const;

        // optimization methods

        bool vertexCoverValidityEdgescheckBitArray(bool* nodeSubset);

        bool vertexCoverValidityNodesRemoved(bool* nodeSubset,const std::vector<uint>& nodesRemoved,const std::vector<uint>& nodesAdded = std::vector<uint>());

        double costFunction(bool* NodeSubset);

        std::vector<uint> getSharedAdjacentNodes(std::vector<uint>& nodes);

        std::vector<uint> getSwappablesIn(std::vector<uint>& nodes,bool* solution);
        std::vector<uint> getSwappablesOut(std::vector<uint>& nodes,bool* solution);

        uint getMaxDegree()const;
        double getAverageDegree()const;

};