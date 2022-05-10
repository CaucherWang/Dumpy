//
// Created by Zeyu Wang on 2021/9/24.
//

#include "../../include/Searchers/DumpyInMemorySearch.h"
#include "../../include/Utils/MathUtil.h"
#include <algorithm>
#include <unordered_set>
#include <random>
#include <chrono>

static float *t_paa;

struct cand_tardis{
    FullAryTreeNode* node;
    double lb;

    cand_tardis(){;}
    cand_tardis(FullAryTreeNode*n, double l){
        node = n;
        lb = l;
    }
};

bool comp_tardis_cand(const cand_tardis& x, const cand_tardis& y){
    return x.lb < y.lb;
}

vector<PqItemSeries *> * DumpyInMemorySearch::approxSearch(FullAryTreeNode* root, float* query, int k, int threshold){
    auto* queryTs = new TimeSeries(query);
    t_paa = queryTs->paa;
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    approxSearchSub(root, sax, queryTs, k, threshold, heap);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

void DumpyInMemorySearch::approxSearchSub(FullAryTreeNode* root, unsigned short *sax, TimeSeries* queryTs, int k, int threshold, vector<PqItemSeries*>*heap){
    if(root->isLeafNode())  return;
    FullAryTreeNode *cur = root->route(sax), *parent = root;
    while (cur!= nullptr && !cur->isLeafNode() && cur->size > threshold) {
        parent = cur;
        cur = cur->route(sax);
    }
    int rest = threshold;

    if(cur != nullptr) {
        cur->exactSearchKnnInMemory(k, queryTs, *heap);
        rest -= cur->size;
    }

    vector<cand_tardis>candidates;
    candidates.reserve((*parent->children).size());
    for(auto &iter:(*parent->children))
        if(iter.second != nullptr && iter.second != cur)
            candidates.emplace_back(iter.second, SaxUtil::LowerBound_Paa_iSax(t_paa, iter.second->sax, iter.second->layer));

    sort(candidates.begin(), candidates.end(), comp_tardis_cand);

    for(int i=0;i<candidates.size() && rest > 0;++i){
        if(rest >= candidates[i].node->size){
            candidates[i].node->exactSearchKnnInMemory(k, queryTs, *heap);
            rest -= candidates[i].node->size;
        }
        else {
            if(candidates[i].node->isLeafNode())
                candidates[i].node->exactSearchKnnInMemory(k, queryTs, *heap, rest);
            else
                approxSearchSub(candidates[i].node, sax, queryTs, k, rest, heap);
            return;
        }
    }
}


