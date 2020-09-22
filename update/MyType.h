#ifndef MYTYPE_H
#define MYTYPE_H

#include "time.h"
#include <assert.h> 

#define INFINITY 999999
#define EPSILON 1.0E-10f
#define MAXINFO 500000

#define ALPHA 10
#define WEIGHT_RANGE 10000
#define WEIGHT_SETTING 1//zipf(ALPHA,WEIGHT_RANGE)
#define ROADNET_NAME "col.roadnet"
//MAXINFO是路网数据中顶点的最大编号

#define Min(a,b) ((a)<(b))?(a):(b)
#define Max(a,b) ((a)>(b))?(a):(b)

typedef int LENGTH;
typedef double WEIGHT;

enum
{
	TYPE_OF_OTHER = 0,
	TYPE_OF_CUSTOMER = 1,
	TYPE_OF_QUERY = 2
};

LENGTH MyRandom(LENGTH max);
//#define TEST_PROB
#define USE_EC// if you want to use effective KNLC
#endif
