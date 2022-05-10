//
// Created by Zeyu Wang on 2021/8/7.
//

#include "../../include/Utils/TimeSeriesUtil.h"
#include "../../include/Utils/FileUtil.h"
#include <bitset>

using namespace std;

void TimeSeriesUtil::heap_data_copy(vector<PqItemSeries *> &heap){
    for(auto *pis:heap)
        pis->copyData();
}

bool TimeSeriesUtil::isSame(const PqItemSeries *ts1, const float *ts2)
{

    float *t1 = (ts1->ts);
    for(int i = 0; i < Const::tsLength; ++i) {
        if(abs(t1[i] - ts2[i])>1e-5)
            return false;
    }
    return true;;
}

int TimeSeriesUtil::intersectionTsSets(const vector<PqItemSeries *> *tsSet1, vector<float *> *tsSet2){
    int intersectionNum = 0, size= tsSet1->size();
    for(const PqItemSeries* currentTs : *tsSet1)
        for(int i=0;i<size;++i){
            if(isSame(currentTs, (*tsSet2)[i])) {
                intersectionNum++;
                break;
            }
        }

    return intersectionNum;
}

int TimeSeriesUtil::intersectionTsSetsCardinality(const vector<PqItemSeries *> &tsSet1, const vector<PqItemSeries *> &tsSet2){
    int intersectionNum = 0;
    for(PqItemSeries* currentTs : tsSet1)
        for (PqItemSeries* targetTs : tsSet2)
            if (isSame(currentTs, targetTs->ts)) {
                intersectionNum += 1;
                break;
            }
    return intersectionNum;
}


double TimeSeriesUtil::euclideanDist(const float* ts_1, const float* ts_2, int len) {
    double sum = 0, dp;
    for (int i = 0; i < len; i++) {
        dp = ts_1[i] - ts_2[i];
        sum += dp * dp;
    }
    return sum;
}

double TimeSeriesUtil::euclideanDist(const float* ts_1, const float* ts_2, int len, double bound) {
    double sum = 0, dp;
    for (int i = 0; i < len && sum < bound; i++) {
        dp = ts_1[i] - ts_2[i];
        sum += dp * dp;
    }
    return sum;
}

template<typename ... Args>
string TimeSeriesUtil::str_format(const string &format, Args ... args)
{
    auto size_buf = std::snprintf(nullptr, 0, format.c_str(), args ...) + 1;
    std::unique_ptr<char[]> buf(new(std::nothrow) char[size_buf]);

    if (!buf)
        return string{};

    std::snprintf(buf.get(), size_buf, format.c_str(), args ...);
    return std::string(buf.get(), buf.get() + size_buf - 1);
}