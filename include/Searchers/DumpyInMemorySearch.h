//
// Created by Zeyu Wang on 2021/9/24.
//

#ifndef DUMPY_DUMPYINMEMORYSEARCH_H
#define DUMPY_DUMPYINMEMORYSEARCH_H


#include <vector>
#include "../DataStructures/PqItemSeries.h"
#include "../DataStructures/FullAryTreeNode.h"

class DumpyInMemorySearch {
public:
    static vector<PqItemSeries *> *approxSearch(FullAryTreeNode *root, float *query, int k, int threshold);

    static void approxSearchSub(FullAryTreeNode *root, unsigned short *sax, TimeSeries *queryTs, int k, int threshold,
                                vector<PqItemSeries *> *heap);

};


#endif //DUMPY_DUMPYINMEMORYSEARCH_H
