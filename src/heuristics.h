
#include "utilities.h"
#include "WeightedVertexGraph.h"
#include <cstddef>
#include <algorithm>
#include <unordered_set>
#include <utility>


bool vertexCoverValidityEdgescheckBit(WeightedVertexGraph* graph,NodeBitList* NodeSubset);

bool vertexCoverValidityEdgescheck(WeightedVertexGraph* graph,NodeList* NodeSubset);

double costFunction(WeightedVertexGraph* graph,NodeList* NodeSubset);
double costFunctionBit(WeightedVertexGraph* graph,NodeBitList* NodeSubset);

NodeBitList* greedySolution(WeightedVertexGraph* graph);