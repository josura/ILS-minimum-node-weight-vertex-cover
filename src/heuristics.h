
#include "utilities.h"
#include "WeightedVertexGraph.h"
#include <cstddef>

//do not work
bool vertexCoverValidityNodescheck(WeightedVertexGraph* graph,NodeList* NodeSubset);
bool vertexCoverValidityNodescheckBit(WeightedVertexGraph* graph,NodeBitList* NodeSubset);


bool vertexCoverValidityEdgescheckBit(WeightedVertexGraph* graph,NodeBitList* NodeSubset);

bool vertexCoverValidityEdgescheck(WeightedVertexGraph* graph,NodeList* NodeSubset);

double costFunction(WeightedVertexGraph* graph,NodeList* NodeSubset);
double costFunctionBit(WeightedVertexGraph* graph,NodeBitList* NodeSubset);