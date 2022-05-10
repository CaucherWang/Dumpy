//
// Created by Zeyu Wang on 2022/1/13.
//

#include <thread>
#include <cassert>
#include <cmath>
#include "../../include/DataStructures/DumpyNode.h"
#include "../../include/Utils/FileUtil.h"
#include "../../include/Utils/MathUtil.h"
#include "../../include/Utils/TimeSeriesUtil.h"
#include "../../include/Utils/SaxUtil.h"

unsigned short *DumpyNode::saxes = nullptr;
float *DumpyNode::paas = nullptr;
int DumpyNode::a = MathUtil::nChooseK(Const::segmentNum, 1), DumpyNode::b = MathUtil::nChooseK(Const::segmentNum, 2), DumpyNode::c = MathUtil::nChooseK(Const::segmentNum, 3);
const int DumpyNode::power_2[]{1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536};
extern long SAX_PAA_CPU_TIME , SAX_PAA_READ_TIME ;
static long MAT1_CPU_TIME_STAT = 0, MAT1_CPU_TIME_COPY = 0, MAT1_WRITE_TIME = 0, MAT1_TOTAL_TIME = 0, MAT1_READ_TIME = 0,
        MAT2_CPU_TIME = 0, MAT2_WRITE_TIME = 0, MAT2_TOTAL_TIME = 0, SAX_PAA_TOTAL_TIME = 0, MAT2_READ_TIME = 0,
        GROW_CPU_TIME = 0,  GROW_CPU_TIME_1st = 0, GROW_TOTAL_TIME = 0, SMALL_FILES_BYTES_WRITE = 0, SMALL_FILES_BYTES_READ = 0,
        RAND_READ_CNT = 0, RAND_WRITE_CNT = 0, SEQ_READ_CNT = 0, SEQ_WRITE_CNT = 0;

DumpyNode* DumpyNode::route(const unsigned short *_sax){
    if(isLeafNode())
        return this;
    int nav_id = SaxUtil::extendSax(_sax, bits_cardinality, chosenSegments);
    if(children[nav_id] == nullptr) return this;
    return children[nav_id]->route(_sax);
}

DumpyNode* DumpyNode::route1step(const unsigned short *_sax){
    assert(!isLeafNode());
    int nav_id;
    if(layer >= 1)
        nav_id = SaxUtil::extendSax(_sax, bits_cardinality, chosenSegments);
    else
        nav_id = SaxUtil::invSaxHeadFromSax(_sax, Const::bitsCardinality, Const::segmentNum);
    return children[nav_id];
}

void DumpyNode::search(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir) const{
    assert(isLeafNode());
    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
    string fn = index_dir+ getFileName();

//    long fs = FileUtil::getFileSize(fn.c_str());
//    int series_num = fs / Const::tsLengthBytes;
//    assert(series_num == size);

    FILE *f = fopen(fn.c_str(), "rb");
    struct timeval io{};
    Const::timer_start(&io);
    auto *ts = new float[size * Const::tsLength];
//    for(int i=0;i<size;++i)
//        fread(ts + i * Const::tsLength, sizeof(float), Const::tsLength, f);
    fread(ts, sizeof(float), size * Const::tsLength, f);

    for(int i=0;i<size;++i){
//        struct timeval start{};
//        Const::timer_start(&start);
        double dist = TimeSeriesUtil::euclideanDist(queryTs->ts, ts + i * Const::tsLength, Const::tsLength, bsf);
//        DIST_CALC_TIME += Const::timer_end(&start);

        if(heap.size() < k){
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }else if(dist < bsf){
            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            delete heap.back();
            heap.pop_back();
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }

        if(heap.size() >= k)    bsf = heap[0]->dist;
    }

    for(PqItemSeries*s: heap){
        if(s->needDeepCopy) s->copyData();
    }
    delete[]ts;
    fclose(f);
}

void DumpyNode::search(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir, unordered_set<float*, createhash, isEqual>*hash_set) const{
    assert(isLeafNode());
    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
    string fn = index_dir+ getFileName();

    FILE *f = fopen(fn.c_str(), "rb");
    auto *ts = new float[size * Const::tsLength];
    fread(ts, sizeof(float), size * Const::tsLength, f);

    for(int i=0;i<size;++i){
        if(hash_set->find(ts + i * Const::tsLength) != hash_set->end()) continue;
        double dist = TimeSeriesUtil::euclideanDist(queryTs->ts, ts + i * Const::tsLength, Const::tsLength, bsf);

        if(heap.size() < k){
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            hash_set->insert(ts + i * Const::tsLength);
        }else if(dist < bsf){
            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            hash_set->erase(heap.back()->ts);
            delete heap.back();
            heap.pop_back();
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            hash_set->insert(ts + i * Const::tsLength);
        }

        if(heap.size() >= k)    bsf = heap[0]->dist;
    }

    for(PqItemSeries*s: heap){
        if(s->needDeepCopy) s->copyData();
    }
    delete[]ts;
    fclose(f);
}

// put actual series into disk file of nodes in 1st layer
void materialize1stLayer(string datafn, DumpyNode* root, int *navids, string index_dir){
    auto start_t = chrono::system_clock::now();
    Const::logPrint("Start move data to disk file in 1st layer.");
    FILE *f = fopen(datafn.c_str(), "r");
    long rest = root->size, total = root->size, cur = 0;
    unordered_map<DumpyNode*, FBL_UNIT>fbl;

    RAND_READ_CNT++;
    SEQ_READ_CNT += rest;

    // There is another implementation method that fbl in each node stores a pointer vector where each pointer points to a series.
    // This may incur many write calls.
    while(rest > 0){
        fbl.clear();
        long num;
        if(rest > Const::fbl_series_num)    num = Const::fbl_series_num;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        auto end = chrono::system_clock::now();
        fread(tss, sizeof(float),num * Const::tsLength,  f);
        auto start = chrono::system_clock::now();
        MAT1_READ_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

        // statistic the size of each node fbl size, and allocate memory for them
        for(long i=cur;i<cur+num;++i)
            fbl[root->children[navids[i]]].size++;
        for(auto &iter:fbl)
            iter.second.buffer = new float [(long)iter.second.size * Const::tsLength];
        end = chrono::system_clock::now();
        MAT1_CPU_TIME_STAT += chrono::duration_cast<chrono::microseconds>(end - start).count();

        //copy series to node to ensure write is only called once upon a node fbl
        for(long i = cur; i<cur+num;++i){
            FBL_UNIT* fbl_node = &fbl[root->children[navids[i]]];
            copy(tss + (i-cur) * Const::tsLength, tss + (i+1-cur)* Const::tsLength, fbl_node->buffer + (long)fbl_node->pos++ * Const::tsLength);
        }
        start = chrono::system_clock::now();
        MAT1_CPU_TIME_COPY += chrono::duration_cast<chrono::microseconds>(start - end).count();

        start = chrono::system_clock::now();
        // write series in order to node file from node fbl
        for(auto & iter:fbl){
            string outfile = index_dir ;
            int id= iter.first->id;
            if(iter.first->partition_id == -1)
                outfile += "U_" + to_string(iter.first->id);
            else
                outfile += to_string(iter.first->layer) + "_" + to_string(iter.first->partition_id);
            FILE *outf = fopen(outfile.c_str(), "a");

            RAND_WRITE_CNT++;
            SEQ_WRITE_CNT += iter.second.size;
//            long bytes = iter.second.size * Const::tsLengthBytes;
//            if(bytes >= Const::small_file_threshold)  SMALL_FILES_BYTES_WRITE += bytes;

            fwrite(iter.second.buffer, sizeof(float), iter.second.size * Const::tsLength, outf);
            fclose(outf);
            delete[]iter.second.buffer;
        }
        end = chrono::system_clock::now();
        MAT1_WRITE_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();
        delete[] tss;

        rest-=num;
        cur += num;
        Const::logPrint("Now in 1st layer " + to_string((double)cur / (double)total * 100) + "% series have been written to disk.");

    }

    fclose(f);
    delete[] navids;
    auto end_t = chrono::system_clock::now();
    MAT1_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(end_t - start_t).count();
}

// put actual series into disk file of nodes below 1st layer from file in 1st layer
void materializeInterNode(DumpyNode *node, unsigned short *saxes) {
    auto start_t = chrono::system_clock::now();

    FILE *f = fopen((Const::idxfn + "U_" + to_string(node->id)).c_str(), "r");
    long rest = node->size, cur = 0, num;
    unordered_map<DumpyNode*, LBL_UNIT>lbl;

//    long bytes = rest * Const::tsLengthBytes;
//    if(bytes >= Const::small_file_threshold)    SMALL_FILES_BYTES_READ += bytes;

    RAND_READ_CNT++;
    SEQ_READ_CNT+=rest;

    while(rest > 0){
        lbl.clear();
        if(rest > Const::fbl_series_num)    num = Const::fbl_series_num;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        auto end = chrono::system_clock::now();
        fread(tss, sizeof(float),num * Const::tsLength,  f);
        auto start = chrono::system_clock::now();
        MAT2_READ_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

        for(long i = cur; i < cur + num; ++i){
            DumpyNode* target = node->route(saxes + (long)node->offsets[i] * Const::segmentNum);
            lbl[target].buffer.push_back(tss + (i - cur) * Const::tsLength);
        }

        end = chrono::system_clock::now();
        MAT2_CPU_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();


        for(auto &iter:lbl){
            string outfile = Const::idxfn + iter.first->getFileName();
            FILE *outf = fopen(outfile.c_str(), "a");
//            SEQ_WRITE_CNT += iter.second.buffer.size();
//            bytes = iter.second.buffer.size() * Const::tsLengthBytes;
//            if(bytes >= Const::small_file_threshold)    SMALL_FILES_BYTES_WRITE += bytes;
            RAND_WRITE_CNT++;
            for(float *dat:iter.second.buffer)
                fwrite(dat, sizeof(float), Const::tsLength, outf);
            fclose(outf);
        }
        start = chrono::system_clock::now();
        MAT2_WRITE_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

        delete[]tss;
        rest-=num;
        cur += num;
    }
    fclose(f);
    FileUtil::FileRemove((Const::idxfn + "U_" + to_string(node->id)).c_str());
    auto end_t = chrono::system_clock::now();
    MAT2_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(end_t - start_t).count();
}

DumpyNode *DumpyNode::BuildIndex(string &datafn, string &saxfn, string &paafn, vector<vector<int>> *g) {
    Const::logPrint("Start building index.");
    auto start_t = chrono::system_clock::now();
    FileUtil::checkDirClean(Const::idxfn.c_str());
    auto end = chrono::system_clock::now();
    long series_num = generateSaxAndPaaTbl();
//    int series_num = loadSax(saxfn);
//    loadPaa(paafn);
    auto start = chrono::system_clock::now();
    Const::logPrint("Finish building sax and paa table.");
    SAX_PAA_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

    auto* root = new DumpyNode();
    root->size = series_num;
    for(int &i:root->bits_cardinality)  i=0;
    partUnit nodeIn1stLayer[Const::vertexNum];
    int *navids = new int[series_num];
    for(int i=0;i<Const::vertexNum;++i)
        nodeIn1stLayer[i].id = i, nodeIn1stLayer[i].size=0, nodeIn1stLayer[i].pid = -1;

    // get 1st layer node size
    for(long i=0;i<series_num;++i){
        unsigned short *asax = saxes + i * Const::segmentNum;
        int nav_id = SaxUtil::invSaxHeadFromSax(asax, Const::bitsCardinality, Const::segmentNum);
        navids[i] = nav_id;
        nodeIn1stLayer[nav_id].size++;
    }

    Const::logPrint("Finish statistic size of nodes in the 1st layer.");

    // partition 1st layer
    int partNum = partition1stLayer(nodeIn1stLayer, g, Const::alpha_1st);
    Const::logPrint("Finish partition");
    DumpyNode* childrenList[partNum];
    for(int i=0;i<partNum;++i)  childrenList[i] = new DumpyNode(1, i);
    root->children.resize(Const::vertexNum);
    for(int i=0;i<Const::vertexNum;++i){
        if(nodeIn1stLayer[i].pid == -1) {
            assert(nodeIn1stLayer[i].size > Const::th);
            root->children[i] = new DumpyNode(1, nodeIn1stLayer[i].size, i);
            root->children[i]->generateSaxAndCardIn1stLayer(i);
        }
        else{
            int pid = nodeIn1stLayer[i].pid;
            root->children[i] = childrenList[pid];
            childrenList[pid]->size += nodeIn1stLayer[i].size;
            childrenList[pid]->generateSaxAndCardIn1stLayer4LeafNode(i);
        }
    }
    Const::logPrint("Finish build index structure 1st layer.");

    thread IO(materialize1stLayer, datafn, root, navids, Const::idxfn);
    
    // put data offsets to internal nodes in 1st layer
    for(int i=0;i<Const::vertexNum;++i)
        if(nodeIn1stLayer[i].size > Const::th)
            root->children[i]->offsets.reserve(nodeIn1stLayer[i].size);
    for(int i=0;i<series_num;++i){
        int nav_id = navids[i];
        if(nodeIn1stLayer[nav_id].size > Const::th) {
            root->children[nav_id]->offsets.push_back(i);
        }
    }
    Const::logPrint("data offsets have been put into nodes in the 1st layer.");
    end = chrono::system_clock::now();
    GROW_CPU_TIME_1st += chrono::duration_cast<chrono::microseconds>(end - start).count();

    int j = 0;
    Const::logPrint("start grow the index structure");
    for(int i=0;i<Const::vertexNum;++i){
        if(nodeIn1stLayer[i].size > Const::th)
            root->children[i]->growIndex();
        if(++j%10000 == 0)  Const::logPrint(to_string(j) + " nodes in the 1st layer has been processed.");
    }
    start = chrono::system_clock::now();
    GROW_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

    Const::logPrint("build index skeleton finished.");

    IO.join();
    Const::logPrint("Start materialize internal nodes in the 1st layer");
    for(int i=0;i<Const::vertexNum;++i)
        if(nodeIn1stLayer[i].size > Const::th)
            materializeInterNode(root->children[i], saxes);
    Const::logPrint("build index successfully!");
    auto end_t = chrono::system_clock::now();
    cout << "Total time is " << chrono::duration_cast<chrono::microseconds>(end_t - start_t).count() / 1000 << "ms."<<endl;
    cout << "Building sax and paa total time is "<<SAX_PAA_TOTAL_TIME / 1000 <<"ms, cpu time is "
    << SAX_PAA_CPU_TIME / 1000 <<"ms, I/O read time is " << SAX_PAA_READ_TIME / 1000 << "ms."<<endl;

    cout << "During the process of building 1st layer index structure, CPU time is "<< GROW_CPU_TIME_1st / 1000 <<"ms, "<<endl;

    cout << "During the process of materializing 1st layer nodes, total time is "<< MAT1_TOTAL_TIME / 1000
         <<"ms, I/O read time is "<< MAT1_READ_TIME / 1000 <<"ms, CPU statistic time  is " << MAT1_CPU_TIME_STAT / 1000
         << "ms, CPU copy Time is " << MAT1_CPU_TIME_COPY / 1000 << "ms, while I/O write time  is " << MAT1_WRITE_TIME / 1000 << "ms. "<< endl;

    cout << "During the process of growing index structure, total time is "<< GROW_TOTAL_TIME / 1000
        <<"ms, other CPU time is "<< GROW_CPU_TIME / 1000 <<"ms."<<endl;

    cout << "During the process of materializing internal nodes, total time is "<<MAT2_TOTAL_TIME / 1000
         <<"ms, I/O read time is "<<MAT2_READ_TIME / 1000 <<"ms, CPU time is " << MAT2_CPU_TIME / 1000
         << "ms, while I/O write time is " << MAT2_WRITE_TIME / 1000 << "ms." << endl;

    cout << "Random read count = " << RAND_READ_CNT << endl
        << "Random write count = " << RAND_WRITE_CNT << endl
        <<"Sequence read count = " << SEQ_READ_CNT << endl
        <<"Sequence write count = " << SEQ_WRITE_CNT << endl
        <<"Small file read bytes = "<<SMALL_FILES_BYTES_READ << endl
        <<"Small file write bytes = " << SMALL_FILES_BYTES_WRITE << endl;

    return root;
}

void DumpyNode::growIndex() {
    if(size <= Const::th)   return;
    auto start = chrono::system_clock::now();
    int chosen_num = SaxUtil::findFirstGE(power_2, 1, Const::segmentNum + 1, size / Const::th + 1);
    PAA_INFO* paa = statPaa();
    chooseSegment(paa, chosen_num);
    
    // statistic children information in order to partition
    partUnit nodes[1<<chosen_num];
    for(int i=0;i<(1<<chosen_num);++i)
        nodes[i].id = i, nodes[i].size=0, nodes[i].pid = -1;
    vector<vector<int>>node_offsets(1<<chosen_num, vector<int>());

    for(int i=0;i<size;++i){
        int new_id = SaxUtil::extendSax(DumpyNode::saxes + (long)offsets[i] * (Const::segmentNum), bits_cardinality, chosenSegments);
        nodes[new_id].size++;
        node_offsets[new_id].push_back(offsets[i]);
    }

    if(this->layer > 1) vector<int>().swap(offsets);

    int partNum = partition(nodes);
    // build rest data node if any
    for(auto &node:nodes)
        if(node.size <= Const::th && node.pid == -1)
            node.pid = ++partNum;
        
    
    DumpyNode* childrenList[partNum];
    for(int i=0;i<partNum;++i)  childrenList[i] = new DumpyNode(this, i);
    children.resize(1 << chosen_num);
    for(int i=0;i<(1 << chosen_num);++i){
        if(nodes[i].size <= 0)  continue;
        else if(nodes[i].pid == -1) {
            children[i] = new DumpyNode(this, nodes[i].size, i);
            generateSaxAndCardinality(children[i], i);
            children[i]->offsets.resize(nodes[i].size);
            copy(node_offsets[i].begin(),  node_offsets[i].end(), children[i]->offsets.begin());
            vector<int>().swap(node_offsets[i]);
        }
        else{
            int _pid = nodes[i].pid;
            children[i] = childrenList[_pid];
            childrenList[_pid]->size += nodes[i].size;
            generateSaxAndCardinality4LeafNode(children[i], i);
            vector<int>().swap(node_offsets[i]);
        }
    }
    
    vector<vector<int>>().swap(node_offsets);
    auto end = chrono::system_clock::now();
    GROW_CPU_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();
    
    for(auto &child: children){
        if(child!= nullptr && child->size > Const::th){
            child->growIndex();
        }
    }

}

// TODO: this function may be optimized with SIMD
PAA_INFO* DumpyNode::statPaa(){
    auto* r = new PAA_INFO();
    double split_line[Const::segmentNum],paa_max[Const::segmentNum],paa_min[Const::segmentNum],paa_mu[Const::segmentNum];
    // TODO: optimize
    double lb;  // ub is the new split line
    for(int i=0; i < Const::segmentNum; ++i)
        SaxUtil::getValueRange(sax[i] << 1, bits_cardinality[i] + 1, &lb, &split_line[i]);
    for(auto &i:paa_max) i = - numeric_limits<double>::max();
    for(auto &i:paa_min) i = numeric_limits<double>::max();
    for(auto &i:r->paa_up_size) i = 0;
    for(auto &i:r->paa_below_size) i = 0;
    for(auto &i:r->paa_variance) i = 0;
    for(auto &i:paa_mu) i=0;
    for(long offset:offsets){
        float * start = paas + offset * (Const::segmentNum);
        for(int i=0;i<Const::segmentNum;++i){
            double value = *(start + i);
            paa_mu[i] += value;
            paa_min[i] = min(paa_min[i], value);
            paa_max[i] = max(paa_max[i], value);
            if(value > split_line[i]) {
                r->paa_up_size[i]++;
            }
            else {
                r->paa_below_size[i]++;
            }
        }
    }
    for(double & i : paa_mu) {
        i /= size;
    }

    for(long offset:offsets){
        float * start = paas + offset * (Const::segmentNum);
        for(int i=0;i<Const::segmentNum;++i){
            double value = *(start + i);
            r->paa_variance[i] += (value - paa_mu[i]) * (value - paa_mu[i]);
        }
    }
    return r;
}

struct tmp{
    int i{};
    double score{};
    tmp(int _i, double _score){i=_i;score = _score;}
    tmp(){;}

    static bool order(tmp a,tmp b){
        return a.score < b.score;
    }

    static bool orderdesc(tmp a,tmp b){
        return a.score > b.score;
    }
};

int DumpyNode::chooseOneSegment(PAA_INFO* node){
    int min = numeric_limits<int>::max(), min_index = -1;
    for(int i = 0; i < Const::segmentNum;++i){
        int big = max(node->paa_up_size[i], node->paa_below_size[i]);
        if(big < min){
            min = big;
            min_index = i;
        }
    }
    return min_index;
}

void DumpyNode::chooseSegment(PAA_INFO *paa, int chosen_num) {
    chosenSegments.resize(chosen_num);
    if(chosen_num == 1) { chosenSegments[0]=chooseOneSegment(paa) ; return;}

    tmp scores[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        if(bits_cardinality[i] >= Const::bitsCardinality)
            scores[i] = tmp(i, -1);
        else
            scores[i] = tmp(i, paa->paa_variance[i]);
    sort(scores, scores+Const::segmentNum, tmp::orderdesc);
    for(int i=0;i<chosen_num;++i)
        chosenSegments[i] = scores[i].i;
    sort(chosenSegments.begin(), chosenSegments.end());
}

void DumpyNode::generateSaxAndCardIn1stLayer(int new_id){
    for(int i = Const::segmentNum - 1; i >=0 ;--i){
        unsigned short t = new_id % 2 ;
        new_id >>= 1;
        sax[i] =  t;
    }
}

void DumpyNode::generateSaxAndCardinality(DumpyNode* node, int new_id){
    copy(sax, sax + Const::segmentNum, node->sax);
    copy(bits_cardinality, bits_cardinality + Const::segmentNum, node->bits_cardinality);
    for(int i = chosenSegments.size() - 1; i >=0 ;--i){
        int seg = (chosenSegments)[i];
        node->bits_cardinality[seg]++;
        int t = new_id % 2 ;
        new_id >>= 1;
        node->sax[seg] = (node->sax[seg] << 1) + t;
    }
}

void DumpyNode::generateSaxAndCardIn1stLayer4LeafNode(int new_id){
    if(bits_cardinality[0] == -1){
        for(int &i:bits_cardinality)    i=1;
        generateSaxAndCardIn1stLayer(new_id);
        return;
    }
    for(int i = Const::segmentNum - 1; i >=0 ;--i){
        int t = new_id % 2 ;
        new_id >>= 1;
        if(bits_cardinality[i] == 1 && sax[i] != t){
            bits_cardinality[i] = 0;
        }
    }
}

void DumpyNode::generateSaxAndCardinality4LeafNode(DumpyNode* node, int new_id){
    if(node->bits_cardinality[0] == -1){
        generateSaxAndCardinality(node, new_id);
        return;
    }
    for(int i = chosenSegments.size() - 1; i >=0 ;--i){
        int seg = chosenSegments[i];
        int t = new_id % 2 ;
        new_id >>= 1;
        if(node->bits_cardinality[seg] == bits_cardinality[seg] + 1 && node->sax[seg] % 2 != t){
            node->bits_cardinality[seg]--;
            node->sax[seg] >>= 1;
        }
    }
}

struct failUnit{
    partUnit *node{};
    int neighbor_size{};

    failUnit(partUnit* a, int b){ node = a; neighbor_size = b;}
};

static bool comp_fail_node(failUnit *x, failUnit *y){
    return x->neighbor_size > y->neighbor_size;
}

int DumpyNode::partition1stLayer(partUnit *nodes_map, vector<vector<int>> *g, double filling_factor) {
    vector<partUnit*>nodes;
    int total_size = 0, node_number;
    for(int i=0;i<Const::vertexNum;++i){
        if(nodes_map[i].size <= Const::th && nodes_map[i].size > 0) {
            nodes.push_back(&nodes_map[i]); 
            total_size += nodes_map[i].size;
        }
    }

    node_number = nodes.size();
    if(node_number < 2) return 0;
    if(node_number == 2 && nodes[0]->size + nodes[1]->size > Const::th) return 0;
//    cout << "node number = " << node_number <<", total size = " << total_size << endl;
    sort(nodes.begin(),  nodes.end(), partUnit::comp_size);
    int pid = 0;
    int k = total_size / Const::th + 1;

    int _id, finish_num = 0, finish_size = 0, temp_size = 0, fail_partition_size = 0;
    partUnit *temp_node;
    vector<partUnit*>temp, candidates;
    // first loop, build 2-clique
    {
        for(int i=0;i<node_number;++i){
            if(nodes[i]->pid != -1 || nodes[i]->size > Const::th)   continue;
            temp_size = nodes[i]->size;  temp.clear(); candidates.clear();  temp.push_back(nodes[i]);
            _id = nodes[i]->id;
            for(int j=0; j<a; ++j)
            {
                temp_node = &nodes_map[(*g)[_id][j]];
                if(temp_node->size > Const::th)    continue;
                if(temp_node->pid == -1 && temp_node->size < Const::th)
                    candidates.push_back(temp_node);
            }

            sort(candidates.begin(),candidates.end(), partUnit::comp_size);
            for(partUnit* node:candidates){
                if(node->size + temp_size > Const::th) continue;
                temp.push_back(node);
                temp_size += node->size;
            }

            // fulfill the partition requirement
            if(temp_size >= filling_factor * Const::th){
                nodes[i]->pid = pid;
                for(partUnit* cur: temp)  cur->pid = pid; 
                finish_num += temp.size();
                finish_size += temp_size;
//                cout << pid << ":" << temp_size << endl;
                ++pid;
            }
        }
    }

    if(finish_num >= node_number)   return pid;
//    cout<< "After first loop, finish nodes number = " << finish_num <<", finish series number = " << finish_size <<endl
//    <<"Built " << pid <<" partitions, with avg size = " << accumulate(partition_size.begin(), partition_size.end(), 0.0) / (double)pid<<endl;

    vector<partUnit*>candidates2;
    // second loop, build 4-clique
    {
        for(int i=0;i<node_number;++i){
            if(nodes[i]->pid != -1 || nodes[i]->size > Const::th)   continue;
            _id = nodes[i]->id;
            temp_size = nodes[i]->size;  temp.clear(); temp.push_back(nodes[i]);  candidates2.clear();
            for(int j=0; j<a; ++j) {
                temp_node = &nodes_map[(*g)[_id][j]];
                if(temp_node->pid == -1 && temp_node->size + temp_size <= Const::th){
                    temp.push_back(temp_node);
                    temp_size += temp_node->size;
                }
            }
            for(int j=a;j<a+b;++j){
                temp_node = &nodes_map[(*g)[_id][j]];
                if(temp_node->pid == -1 && temp_node->size < Const::th){
                    candidates2.push_back(temp_node);
                }
            }

            sort(candidates2.begin(),candidates2.end(), partUnit::comp_size);
            for(partUnit* node:candidates2){
                if(node->size + temp_size > Const::th) continue;
                temp.push_back(node);
                temp_size += node->size;
            }

            // fulfill the partition requirement
            if(temp_size >= filling_factor * Const::th){
                nodes[i]->pid = pid;
                for(partUnit* cur: temp)  cur->pid = pid;
                finish_num += temp.size();
                finish_size += temp_size;
                ++pid;
            }
        }
    }

    if(finish_num >= node_number)   return pid;
//    cout<< "After second loop, finish nodes number = " << finish_num <<", finish series number = " << finish_size <<endl
//    <<"Built " << pid <<" partitions, with avg size = " << accumulate(partition_size.begin() + one_loop_pid, partition_size.end(), 0.0) / (double)(pid-one_loop_pid)<<endl;
    vector<failUnit*>fail_partition_node;

    {
        int no_part;
        for(int i=0;i<node_number;++i){
            if(nodes[i]->pid != -1)   continue;
            no_part = 0;
            _id = nodes[i]->id;
            for(int j=0;j<a+b+c;++j){
                temp_node = &nodes_map[(*g)[_id][j]];
                if(temp_node->size > Const::th)    continue;
                if(temp_node->pid == -1){
                    no_part+= temp_node->size;
                }
            }

            fail_partition_node.push_back(new failUnit(nodes[i], no_part + nodes[i]->size));
            fail_partition_size += nodes[i]->size;
        }
    }
//    cout<< "After filling stage, finish nodes number = " << finish_num <<", finish series number = " << finish_size <<endl
//    <<"Built " << pid <<" partitions, with avg size = " << accumulate(partition_size.begin(), partition_size.end(), 0.0) / (double)pid<<endl;
//    cout <<"Fail partition node number = " << fail_partition_node.size() << " , total size = " << fail_partition_size << endl;

    sort(fail_partition_node.begin(),  fail_partition_node.end(), comp_fail_node);
    partUnit* debug_node = &nodes_map[33];
    vector<bool>flag(Const::vertexNum, false);
    for(failUnit* node:fail_partition_node){
        if(flag[node->node->id])    continue;
        temp_size = node->node->size;   temp.clear();   temp.push_back(node->node);
        node->node->pid = pid;  flag[node->node->id] = true;
        finish_num++; finish_size+=node->node->size;
        _id = node->node->id;

        int j;
        for(j=0; j<a +b +c && temp_size < Const::th; ++j)
        {
            temp_node = &nodes_map[(*g)[_id][j]];
            if(temp_node->pid == -1&& temp_node->size + temp_size <= Const::th)
                temp_node->pid = pid, flag[temp_node->id] = true, finish_num++, finish_size += temp_node->size, temp_size+=temp_node->size, temp.push_back(temp_node);
        }

        if(temp_size >= filling_factor * Const::th) {
            ++pid; continue;
        }
        for(int i=4; i <= Const::r_max && temp_size < Const::th; ++i){
            int mask = (1U << i) - 1, r, last1;
            do{
                {
                    int cur = mask ^ node->node->id;
                    temp_node = &nodes_map[cur];
                    if(temp_node->pid == -1 && temp_node->size + temp_size < Const::th)
                        temp_node->pid = pid, flag[temp_node->id] = true, finish_num++, finish_size += temp_node->size, temp_size+=temp_node->size, temp.push_back(temp_node);
                    if(temp_size >= Const::th) break;
                }
                last1 = mask & -mask;
                r = mask + last1;
                mask = (((r ^ mask) >> 2) / last1) | r;
            } while (r < (1 << Const::segmentNum));
        }
        ++pid;
    }

//    cout<< "After stretch stage, finish nodes number = " << finish_num <<", finish series number = " << finish_size <<endl
//    <<"Built " << pid <<" partitions, with avg size = " << accumulate(partition_size.begin(), partition_size.end(), 0.0) / (double)pid<<endl;

    for(auto* _:fail_partition_node)    delete _;
    return pid;
}

int findFirstLE(vector<partUnit*>&nodes, int start, int end, int target){
    int mid, nn = end + 1;
    while (start < end){
        mid = (start + end) / 2;
        if(nodes[mid]->size > target){
            start = mid + 1;
        }else{
            end = mid;
        }
    }
    if(start == end && nodes[start]->size <= target)    return  start;
    return  nn;
}

int DumpyNode::partition(partUnit* nodes_map){
    vector<partUnit*>nodes;
    int total_size = 0, node_number;
    for(int i=0;i<Const::vertexNum;++i){
        if(nodes_map[i].size <= Const::th) {
            nodes.push_back(&nodes_map[i]);
            total_size += nodes_map[i].size;
        }
    }

    node_number = nodes.size();
    if(node_number < 2) return 0;
    if(node_number == 2 && nodes[0]->size + nodes[1]->size > Const::th) return 0;
//    cout << "node number = " << node_number <<", total size = " << total_size << endl;
    sort(nodes.begin(),  nodes.end(), partUnit::comp_size);
    int pid = 0;
    for(int i=0;i<node_number;++i){
        if(nodes[i]->pid != -1) continue;
        int target = Const::th - nodes[i]->size;
        int start = findFirstLE(nodes, i + 1, node_number - 1, target);
        nodes[i]->pid = pid;
        for(int j=start; j < node_number; ++j){
            if(nodes[j]->pid != -1 || nodes[j]->size > target) continue;
            nodes[j]->pid = pid;
            target -= nodes[j]->size;
        }
        ++pid;
    }
    return pid;
}

void DumpyNode::getIndexStats(){
    int total_leaf_node_num = getLeafNodeNum();
    int total_size = getTotalSize();
    cout << "Total size = " << total_size << endl;
    cout <<"Total nodes number = " << getNodeNum() << endl;
    cout << "Leaf node number = " << total_leaf_node_num << endl;
    cout << "1st layer node number = " << get1stLayerNodesNo() <<endl;
    cout << "1st layer internal node number = " << get1stLayerInterNodesNo() << endl;
    cout << "1st layer internal series number = " << get1stLayerInterNodeSeriesNo() << endl;
    cout << "Max. height = " << getMaxHeight() - 1 <<endl;
    cout << "Avg. Height = " << getSumHeight() / (double) total_leaf_node_num << endl;
    cout <<"Avg. Filling Factor = "<< total_size / (double)total_leaf_node_num / Const::th << endl;
    cout << "Bias leaf node ratio = " << (double)getBiasLeafNodeNum() / total_leaf_node_num << endl;
}

int DumpyNode::get1stLayerInterNodesNo(){
    unordered_set<DumpyNode*>node;
    for(DumpyNode* child:children){
        if(child== nullptr || child->size <= Const::th || node.find(child) != node.end())
            continue;
        node.insert(child);
    }
    return node.size();
}

int DumpyNode::get1stLayerNodesNo(){
    unordered_set<DumpyNode*>node;
    for(DumpyNode* child:children){
        if(child== nullptr || node.find(child) != node.end())
            continue;
        node.insert(child);
    }
    return node.size();
}

int DumpyNode::get1stLayerInterNodeSeriesNo(){
    unordered_set<DumpyNode*>node;
    int ret = 0;
    for(DumpyNode* child:children){
        if(child== nullptr || child->size <= Const::th || node.find(child) != node.end())
            continue;
        node.insert(child);
        ret += child->size;
    }
    return ret;
}

int DumpyNode::getMaxHeight(){
    if(isLeafNode())    return 1;
    int max_height = 0;
    unordered_set<DumpyNode*>hash_map;
    for(DumpyNode*child:children){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        max_height = max(child->getMaxHeight(), max_height);
        hash_map.insert(child);
    }
    return max_height + 1;
}

int DumpyNode::getLeafNodeNum(){
    if(isLeafNode())    return 1;
    int sum = 0;
    unordered_set<DumpyNode*>hash_map;
    for(DumpyNode*child:children){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        sum += child->getLeafNodeNum();
        hash_map.insert(child);
    }
    return sum;
}

// isax bias leaf nodes number
int DumpyNode::getBiasLeafNodeNum(){
    if(isLeafNode()){
        int max_b = 0, min_b=8;
        for(int &bc:bits_cardinality){
            max_b = max(max_b, bc);
            min_b = min(min_b, bc);
        }
        return (max_b - min_b >= 4);
    }
    int sum = 0;
    unordered_set<DumpyNode*>hash_map;
    for(DumpyNode*child:children){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        sum += child->getBiasLeafNodeNum();
        hash_map.insert(child);
    }
    return sum;
}

int DumpyNode::getTotalSize(){
    if(isLeafNode())    return size;
    int sum = 0;
    unordered_set<DumpyNode*>hash_map;
    for(DumpyNode*child:children){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        sum += child->getTotalSize();
        hash_map.insert(child);
    }
    return sum;
}

int DumpyNode::getNodeNum(){
    if(isLeafNode())    return 1;
    int sum = 0;
    unordered_set<DumpyNode*>hash_map;
    for(DumpyNode*child:children){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        sum += child->getNodeNum();
        hash_map.insert(child);
    }
    return sum + 1;
}

int DumpyNode::getSumHeight(){
    if(isLeafNode())    return layer;
    int sum_height = 0;
    unordered_set<DumpyNode*>hash_map;
    for(DumpyNode*child:children){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        sum_height += child->getSumHeight();
        hash_map.insert(child);
    }
    return sum_height;
}

int DumpyNode::loadSax(const string & saxfn){
    long f_size = FileUtil::getFileSize(saxfn.c_str()), series_num = f_size / (sizeof(unsigned short) * Const::segmentNum);
    saxes = new unsigned short [f_size / sizeof(unsigned short )];
    FILE *f = fopen(saxfn.c_str(), "rb");
    fread(saxes, sizeof(unsigned short ), f_size / sizeof(unsigned short), f);
    fclose(f);
    Const::logPrint("Finish loading sax");
    return series_num;
}

void DumpyNode::loadPaa(const string & paafn){
    long f_size = FileUtil::getFileSize(paafn.c_str());
    paas = new float [f_size / sizeof(float )];
    FILE *f = fopen(paafn.c_str(), "rb");
    fread(paas, sizeof(float ), f_size / sizeof(float ), f);
    fclose(f);
    Const::logPrint( "Finish loading paa");
}

long DumpyNode::generateSaxAndPaaTbl(){
    string fn = Const::datafn;
    long fs = FileUtil::getFileSize(fn.c_str());
    long series_num = fs / Const::tsLengthBytes;
    cout << "Total Series Number is "<<series_num <<endl;
    float ts[Const::tsLength];
    saxes = new unsigned short[series_num * Const::segmentNum];
    paas = new float[series_num * Const::segmentNum];
    long rest = series_num, cur = 0;
    FILE *f = fopen(fn.c_str(), "rb");

    RAND_READ_CNT++;
    SEQ_READ_CNT += series_num;
    RAND_WRITE_CNT+=2;
    SEQ_WRITE_CNT += series_num;

    while(rest > 0){
        unsigned num;
        if(rest > 4000000)    num = 4000000;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        auto start = chrono::system_clock::now();
        fread(tss, sizeof(float),num * Const::tsLength,  f);
        auto end = chrono::system_clock::now();
        SAX_PAA_READ_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();

        for(int i=0;i<num;++i){
            if(isnan(tss[i * Const::tsLength])){
                for(int j = 0; j < Const::segmentNum; ++j)
                    saxes[i * Const::segmentNum + j] = 0;
                cout << "Dirty data: "<<i << "," <<endl;
            }
            else{
                SaxUtil::paaAndSaxFromTs(tss + i * Const::tsLength,
                                         paas + (cur+ i) * Const::segmentNum, saxes + (cur+ i) * Const::segmentNum,
                                         Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
            }
        }
        delete[] tss;
        rest -=num;
        cur+=num;
        start = chrono::system_clock::now();
        SAX_PAA_CPU_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();
    }

    fclose(f);
    return series_num;
}

DumpyNode *DumpyNode::loadFromDisk(const string &saxfn, const string &idxfn, bool need_sax) {
    if(need_sax)
        loadSax(saxfn);
    ifstream ifs(idxfn, ios::binary);
    boost::archive::binary_iarchive ia(ifs);
    auto *g = new DumpyNode();
    ia >> (*g);
    ifs.close();
    return g;
}
