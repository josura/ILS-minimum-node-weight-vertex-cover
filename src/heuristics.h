
#include "utilities.h"
#include "WeightedVertexGraph.h"
#include <cstddef>
#include <algorithm>
#include <unordered_set>
#include <utility>


bool vertexCoverValidityEdgescheckBitList(WeightedVertexGraph* graph,NodeBitList* NodeSubset);
bool vertexCoverValidityEdgescheckBitArray(WeightedVertexGraph* graph,NodeBitArray NodeSubset);

bool vertexCoverValidityEdgescheck(WeightedVertexGraph* graph,NodeList* NodeSubset);

double costFunction(WeightedVertexGraph* graph,NodeList* NodeSubset);
double costFunctionBitList(WeightedVertexGraph* graph,NodeBitList* NodeSubset);
double costFunctionBitArray(WeightedVertexGraph* graph,NodeBitArray NodeSubset);

NodeBitList* greedySolutionBitList(WeightedVertexGraph* graph);
NodeBitArray greedySolutionBitArray(WeightedVertexGraph* graph);