
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
double costFunction(WeightedVertexGraph* graph,NodeSet* NodeSubset);
double costFunction(WeightedVertexGraph* graph,NodeBitList* NodeSubset);
double costFunction(WeightedVertexGraph* graph,NodeBitArray NodeSubset);

NodeBitList* greedySolutionBitList(WeightedVertexGraph* graph);
NodeSet* greedySolutionBitArray(WeightedVertexGraph* graph,NodeBitArray solution);