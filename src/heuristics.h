
#include "utilities.h"
#include "WeightedVertexGraph.h"
#include <cstddef>


bool vertexCoverValidityEdgescheckBit(WeightedVertexGraph* graph,NodeBitList* NodeSubset);

bool vertexCoverValidityEdgescheck(WeightedVertexGraph* graph,NodeList* NodeSubset);

double costFunction(WeightedVertexGraph* graph,NodeList* NodeSubset);
double costFunctionBit(WeightedVertexGraph* graph,NodeBitList* NodeSubset);