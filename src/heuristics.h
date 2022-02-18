
#include "utilities.h"
#include "WeightedVertexGraph.h"
#include <cstddef>
#include <algorithm>
#include <unordered_set>
#include <utility>
#include <boost/functional/hash.hpp>


bool vertexCoverValidityEdgescheckBitList(WeightedVertexGraph* graph,NodeBitList* NodeSubset);
bool vertexCoverValidityEdgescheckBitArray(WeightedVertexGraph* graph,NodeBitArray NodeSubset);
bool vertexCoverValidityEdgescheckNodeList(WeightedVertexGraph *graph, NodeList& nodeSubset);

bool vertexCoverValidityEdgescheck(WeightedVertexGraph* graph,NodeList* NodeSubset);

double costFunction(WeightedVertexGraph* graph,NodeList& NodeSubset);
double costFunction(WeightedVertexGraph* graph,NodeSet* NodeSubset);
double costFunction(WeightedVertexGraph* graph,NodeBitList* NodeSubset);
double costFunction(WeightedVertexGraph* graph,NodeBitArray NodeSubset);

NodeBitList* greedySolutionBitList(WeightedVertexGraph* graph);
NodeSet* greedySolutionBitArray(WeightedVertexGraph* graph,NodeBitArray solution);

//OPTIMIZATION
NodeBitArray removeRedundantNodes(WeightedVertexGraph *graph, NodeBitArray solution);
NodeBitArray removeRedundantNodesBest(WeightedVertexGraph *graph, NodeBitArray solution);
NodeBitArray removeRedundantNodesSmart(WeightedVertexGraph *graph, NodeBitArray solution);
NodeBitArray randomGreedySolutionAndRemoves(WeightedVertexGraph *graph, NodeBitArray solution);
NodeBitArray randomSmartSolutionAndRemoves(WeightedVertexGraph *graph, NodeBitArray solution);
//OPTIMIZATION

NodeSet* randomGreedySolutionBitArray(WeightedVertexGraph *graph, NodeBitArray solution);
NodeBitArray randomGreedySolutionBitArrayFromPartial(WeightedVertexGraph *graph, NodeBitArray solution);
NodeBitArray randomSmartSolutionBitArrayFromPartial(WeightedVertexGraph *graph, NodeBitArray solution);