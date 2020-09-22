#include<cstdio>
#include<cstring>
#include<iostream>
#include<cstdlib>
#include <sys/time.h>
#include<vector>
#include <xmmintrin.h>
#include<cmath>
#include<map>
#include<cassert>
#include<cstdint>
#include<algorithm>
#include "bitvector.h"
#include "macros.h"

using namespace std;
#define USE_FECTH
#define NO_PRUNE
//#define USE_BIT
//#define USE_BYTE
//#define DX3
//#define REPEAT 100
typedef int LENGTH;
enum NODE_TYPE{
	TYPE_Cut = 3,
	//TYPE_D2 = 2,//leaf
	TYPE_D1 = 1,//single child
	TYPE_C2 = 0,//lca
};
struct timer
{
	timeval stime, etime;

	timer()
	{
		gettimeofday(&stime, NULL);
	}
	void restart()
	{
		gettimeofday(&stime, NULL);
	}

	double getTime()
	{
		gettimeofday(&etime, NULL);
		return (double)(etime.tv_sec - stime.tv_sec) * 1000000 +
			(double)(etime.tv_usec - stime.tv_usec);
	}
};

const LENGTH infinity = 999999999;

int bheight_control;
int *LOG2, *LOGD, **RMQIndex;
int *height;

int topk;
map<int, int> mGoodBranch;
int * vGoodBranch;
double * vProbBranch, * vBetterBranch;
int** goodBranchPos0;

double outputThreshold;
map<int, double> mBidToBetter;

int hash_table[16];

int g_mw;

class NODE
{
public:
	NODE()
	{
		//pos0=pos1=pos2=pos3=pos4=pos5=NULL;
		for(int i = 0; i < 6; i++){
		    positions[i] = nullptr;
		    dlist2[i] = nullptr;
		    dlist3[i] = nullptr;
		    #ifdef USE_BIT
		    bvs[i] = nullptr;
		    #endif
		    #ifdef USE_BYTE
		    bytepos[i]=nullptr;
		    byteSize[i]=0;
		    #endif
		    posSizes[i] = 0;
		}
		//pos0Size=pos1Size=pos2Size=pos3Size=pos4Size=pos5Size=0;
	}

	//int  pos0Size, pos1Size, pos2Size, pos3Size, pos4Size, pos5Size;
	int posSizes[6];
	//int *pos0, *pos1, *pos2, *pos3, *pos4, *pos5;
	int* positions[6];
	LENGTH * dlist2[6];
	LENGTH * dlist3[6];

	int toRMQ;
	#ifdef USE_BIT
	BitVector *bvs[6];
	#endif
	#ifdef USE_BYTE
	uint8_t *bytepos[6];
	int byteSize[6];
	#endif
	~NODE()
	{
		for(int i = 0; i < 6; i++){
		    if(positions[i])
		    {
		        free(positions[i]);
		        positions[i]=nullptr;
		    }
	        if(dlist2[i])
	        {
		        free(dlist2[i]);
		        dlist2[i] = nullptr;
	        }
	        if(dlist3[i])
	        {
	        	free(dlist3[i]);
	        	dlist3[i] = nullptr;
	        }
		    #ifdef USE_BIT
		    if(bvs[i])
		    {
			    delete bvs[i];
			    bvs[i]=nullptr;
		    }
		    #endif
		    #ifdef USE_BYTE
		    if(bytepos[i])
		    {
		    	free(bytepos[i]);
		    	bytepos[i]=nullptr;
		    }
		    #endif
		}
	}

	// Re-initialize the node
	void reinitialize(NODE *vtree, int* goodBranches){

	    // read mToAncestor to position
	    for(int i = 0; i < 5; i++){

	        int mid = i + 1;
	        int branch_id = vGoodBranch[i];

	        const int* lcaPos0 = vtree[branch_id].positions[0];
	        const int lcaPosSize = vtree[branch_id].posSizes[0];

	        int* pos = positions[mid];
	        const int  posSize = posSizes[mid];


	        #ifdef USE_BIT
            if(pos != nullptr && posSize != 0){

                bvs[mid] = new BitVector(lcaPosSize);

                int kk = 0;
                for(int jj = 0; jj < posSize; jj ++){
            //        while(pos[jj] != lcaPos0[kk] && kk < lcaPosSize){
                    while(pos[jj] != lcaPos0[kk]) {
                        kk ++;
                    }

                    bvs[mid]->SetBit(kk);
                    kk++;
                }
            }
            #endif
            #ifdef USE_BYTE
            if(pos != nullptr && posSize != 0){
            	vector<int> vRefIndex(posSize);
            	int index = 0;
            	for(int jj = 0; jj < lcaPosSize && index < posSize; ++jj)
            	{
            		if(lcaPos0[jj] == pos[index])
            		{
            			vRefIndex[index] = jj;
            			++index;
            		}
            	}

                vector<uint8_t> vctemp;
                int last = 0;
                for(int jj = 0; jj < posSize; jj++){
                    if((vRefIndex[jj] - last) < 256)
                    {
                        uint8_t tmp = vRefIndex[jj] - last;
                        if(vRefIndex[jj]-last<0)
                        	assert(false);
                        last = vRefIndex[jj];
                        vctemp.push_back(tmp);
                    }
                    else
                    {
                    	uint8_t tmp = 255;
                    	vctemp.push_back(tmp);
                    	last = last + 255;
                    	--jj;
                        //exit(11);
                    }
                }
                
                byteSize[mid] = vctemp.size();
                bytepos[mid] = (uint8_t*)malloc(sizeof(uint8_t) * byteSize[mid]);
                for(int jj=0; jj<byteSize[mid]; ++jj)
                	bytepos[mid][jj] = vctemp[jj];
            }
            #endif
	    }
	}



	void set(int mid, int size, FILE * fin)
	{
        if(mid >= 0 && mid <= 6){
		    posSizes[mid] = size;
		    positions[mid] = (int*) malloc(sizeof(int) * size);
		    fread(positions[mid], sizeof(int), size, fin);

		} else {
		    int* temp = (int*)malloc(sizeof(int) * size);
		    fread(temp, sizeof(int), size, fin);
		    free(temp);
		}
	}

	void setdlist(NODE * vtree)
	{
		for(int mid=1; mid<6; ++mid)
		{
		    	int bid = vGoodBranch[mid-1];
		    	int * pos = vtree[bid].positions[0];
		    	int psize = vtree[bid].posSizes[0];
		    	if(psize<1)
		    		continue;
		    	int possize = posSizes[mid];
		    	if(possize<1)
		    		continue;
		    	dlist2[mid] = (LENGTH *)malloc(sizeof(LENGTH)* psize);
		    	#ifdef DX3
		    	dlist3[mid] = (LENGTH *)malloc(sizeof(LENGTH)* possize);
		    	assert(positions[mid]);
		    	for(int ii = 0; ii < possize; ++ ii)
		    		dlist3[mid][ii]=dlist2[0][positions[mid][ii]];
		    	#endif
		    	for(int ii = 0; ii < psize; ++ii)
		    	{
		    		dlist2[mid][ii]=dlist2[0][pos[ii]];

		    	}
		    	if(positions[mid])
		    	{
		    		free(positions[mid]);
		    		positions[mid] = nullptr;
		    		#ifdef USE_BYTE
		    		posSizes[mid] = byteSize[mid];
		    		#endif
		    	}
		 }
	}

};

NODE * vtree;

int g_avgw, iVNum;
int iRow;


//int *sizeBranch;
void readANode(FILE * fin)//refer to Node.printMe
{
	//read the vertex id, node type
	int vid = 0, nodetype = 0;
	fread(&vid, sizeof(int), 1, fin);

	NODE & tn = vtree[vid];

	int iBSize = 0;
	fread(&iBSize, sizeof(int), 1, fin);
	if(iBSize > 0)
	{
		//sizePos[vid] = (int *)malloc(sizeof(int)*iBSize);
		//mPos[vid] = (int **)malloc(sizeof(int *)*iBSize);
		//sizeBranch[vid] = iBSize;

		for (int ii = 0; ii < iBSize; ++ii)
		{
			int bid;
			fread(&bid, sizeof(int), 1, fin);
			int mid = mGoodBranch[bid];

			int size = 0;
			fread(&size, sizeof(int), 1, fin);

			if (size>0)
				tn.set(mid, size, fin);
		}
	}

	fread(&nodetype, sizeof(int), 1, fin);

	bool bFake = false;
	int tn_height;
	switch (nodetype)
	{
	case TYPE_D1://single child
	{
		fread(&tn_height, sizeof(int), 1, fin);
	}
	break;
	case TYPE_C2://lca
	{
		//tn.iReal = vid;
		//read the height, used in the pos to do query if the two query vertices have ancestor relationship
		//fread(&tn_branch_height, sizeof(int), 1, fin);
		fread(&tn_height, sizeof(int), 1, fin);

		int iPosSize = 0;
		fread(&iPosSize, sizeof(int), 1, fin);
		tn.set(0, iPosSize, fin);
	}
	break;
	default:
		assert(false);
		break;
	}
	if (!bFake)
	{
		int iDLSize = 0;
		fread(&iDLSize, sizeof(int), 1, fin);

		tn.dlist2[0] = (LENGTH *)malloc(sizeof(LENGTH)*iDLSize);
		fread(tn.dlist2[0], sizeof(LENGTH), iDLSize, fin);
	}

	height[vid] = tn_height;
}

void readBranchInformation(FILE * fin)
{
	fread(&outputThreshold, sizeof(double), 1, fin);//tell which branch we have store prune information
	//read good branches
	fread(&topk, sizeof(int), 1, fin);
	vGoodBranch = (int *)malloc(sizeof(int)*topk);
	fread(vGoodBranch, sizeof(int), topk, fin);//the id, sort by prob in ascending order
	vProbBranch = (double *)malloc(sizeof(double)*topk);
	fread(vProbBranch, sizeof(double), topk, fin);//the prob
	vBetterBranch = (double *)malloc(sizeof(double)*topk);
	fread(vBetterBranch, sizeof(double), topk, fin);//the discount on width

	goodBranchPos0 = (int**) malloc(sizeof(int*) * topk);

	int jj = 0;
	for(int ii = 0; ii < topk; ++ii)
	{
		//cout<<vGoodBranch[ii]<<'\t'<<vProbBranch[ii]<<'\t'<<vBetterBranch[ii]<<endl;
		mBidToBetter[vGoodBranch[ii]] = vBetterBranch[ii];
		//cout << vGoodBranch[ii] << endl;
		hash_table[vGoodBranch[ii] % 14] = ii + 1;
		//if(vBetterBranch[ii] < acceptThreshold && vBetterBranch[ii] < outputThreshold)
		//{
			vGoodBranch[jj] = vGoodBranch[ii];
			++jj;
			mGoodBranch[vGoodBranch[ii]] = jj;


		//}
		//else
		//	mGoodBranch[vGoodBranch[ii]] = -1;
	}
	topk = jj;
	//vGoodBranch.resize(topk);
	//cout<<endl;
	//free(vGoodBranch);
}

void readIndex(string fileIndex)
{
	g_avgw = 0;
	FILE * fin;
	fin = fopen(fileIndex.c_str(), "rb");
	//fopen_s(&fin, fileIndex.c_str(), "rb");

	//read the number of vertices
	fread(&iVNum, sizeof(int), 1, fin);

	vtree = (NODE *)malloc(sizeof(NODE)*iVNum);

	readBranchInformation(fin);


	height = (int *)malloc(sizeof(int)*iVNum);
	//sizeBranch = (int *)malloc(sizeof(int)*iVNum);

	//read all nodes
	for (int ii = 0; ii < iVNum; ++ii)
	{
		readANode(fin);
	}

	for(int ii = 0; ii < iVNum; ++ii){
	    vtree[ii].reinitialize(vtree, vGoodBranch);
	}

	for(int j = 0; j < topk; j++){
	    goodBranchPos0[j] = vtree[vGoodBranch[j]].positions[0];
	}

	//read the part for lca
	int * vtoRMQ = (int *)malloc(sizeof(int)*iVNum);
	fread(vtoRMQ, sizeof(int), iVNum, fin);
	for(int ii = 0; ii < iVNum; ++ii)
		vtree[ii].toRMQ = vtoRMQ[ii];
	free(vtoRMQ);


	fread(&iRow, sizeof(int), 1, fin);
	RMQIndex = (int **)malloc(sizeof(int *) * iRow);
	for (int ii = 0; ii < iRow; ++ii)
	{
		int iCol = 0;
		fread(&iCol, sizeof(int), 1, fin);
		RMQIndex[ii] = (int *)malloc(sizeof(int) *iCol);
		fread(RMQIndex[ii], sizeof(int), iCol, fin);
	}

	for(int ii = 0; ii < iVNum; ++ii)
		vtree[ii].setdlist(vtree);
}

void freeall()
{
	free(height);
	free(vtree);

	free(LOGD);
	free(LOG2);
	for (int ii = 0; ii < iRow; ++ii)
	{
		free(RMQIndex[ii]);
	}
	free(RMQIndex);

	free(vGoodBranch);
}


inline	int LCAQuery(int p, int q){
	//int p = toRMQ[_p], q = toRMQ[_q];

	if (p > q){
		int x = p;
		p = q;
		q = x;
	}
	int len = q - p + 1;

	int i = LOGD[len], k = LOG2[len];

	q = q - i + 1;
	if (height[RMQIndex[k][p]] < height[RMQIndex[k][q]])
		return RMQIndex[k][p];
	else return RMQIndex[k][q];
}

inline	LENGTH H2H(int x, int y){
	if (x == y) return 0;

	NODE & tnx = vtree[x];
	NODE & tny = vtree[y];
	int lca = LCAQuery(tnx.toRMQ, tny.toRMQ);

//	cout << "lca: " << lca << '\t';
	if(lca == x)
	{
	    //++g_avgw;
		return tny.dlist2[0][height[x]];
	}
	else if(lca == y)
	{
		//++g_avgw;
		return tnx.dlist2[0][height[y]];
	}
	else {

		NODE & tnc = vtree[lca];
	//	int hopsize = tnc.pos0Size;
	    int hopsize = tnc.posSizes[0];
	//	int *hop = tnc.pos0;
	    int *hop = tnc.positions[0];

		//g_avgw += hopsize;

		LENGTH *dx = tnx.dlist2[0], *dy = tny.dlist2[0];

		#ifdef USE_FECTH
		_mm_prefetch(dx, _MM_HINT_T0);
		_mm_prefetch(dy, _MM_HINT_T0);
		_mm_prefetch(hop, _MM_HINT_T0);
		#endif

		LENGTH res = infinity;

		//if(hopsize > g_mw)
		//	g_mw = hopsize;
            for (int i = 0; i < hopsize; i++){
            	int j = hop[i];
            	
            	if(j >= height[x])
            		cout<< j;
            	if(j >= height[y])
            		cout<<j;
            	
                LENGTH tmp = dx[j] + dy[j];
                //cout<<j<<'\t'<<tmp<<endl;
                if (res > tmp)
                    res = tmp;
            }
		//cout << res << endl;
		return res;
	}
}


inline	LENGTH P2H(int x, int y){
	if (x == y) return 0;
	const NODE & tnx = vtree[x];
	const NODE & tny = vtree[y];
	const int lca = LCAQuery(tnx.toRMQ, tny.toRMQ);

	if(lca == x)
	{
		//++g_avgw;
		return tny.dlist2[0][height[x]];
	}
	else if(lca == y)
	{
		//++g_avgw;
		return tnx.dlist2[0][height[y]];
	}
	else {

	 //   cout << "lca: " << lca << " ";
		int *hop = nullptr;
		#ifdef USE_BIT
		BitVector* bv = nullptr;
		#endif
		#ifdef USE_BYTE
		uint8_t* bytes = nullptr;
        #endif
		bool bGB=false;
 
        LENGTH *dx=nullptr, *dy=nullptr;
        int hopsize = 0;
        
        for(int mid = 1; mid < 6; ++mid)
        	if(lca == vGoodBranch[mid-1])
        	{
        		#ifdef NO_PRUNE
        		dx = tnx.dlist2[mid];
        		dy= tny.dlist2[mid];
        		#endif
        		#ifdef USE_BIT
        		if(tnx.posSizes[mid] < tny.posSizes[mid])
        		{
	        		bv = tnx.bvs[mid];
	        		hopsize = tnx.posSizes[mid];
	        		#ifdef DX3
	        		dx = tnx.dlist3[mid];
	        		#else
	        		dx = tnx.dlist2[mid];
	        		#endif
	        		dy = tny.dlist2[mid];
        		}
	        	else
	        	{
	        		bv = tny.bvs[mid];
	        		hopsize = tny.posSizes[mid];
	        		#ifdef DX3
	        		dx = tny.dlist3[mid];
	        		#else
	        		dx = tny.dlist2[mid];
	        		#endif
	        		dy = tnx.dlist2[mid];
	        	}
	        	#endif
	        	#ifdef USE_BYTE
	        	if(tnx.posSizes[mid] < tny.posSizes[mid])
	        	{
	        		bytes = tnx.bytepos[mid];
	        		hopsize = tnx.posSizes[mid];
	        		#ifdef DX3
	        		dx = tnx.dlist3[mid];
	        		#else
	        		dx = tnx.dlist2[mid];
	        		#endif
	        		dy = tny.dlist2[mid];
	        	}
	        	else
	        	{
	        		bytes = tny.bytepos[mid];
	        		hopsize = tny.posSizes[mid];
	        		#ifdef DX3
	        		dx = tny.dlist3[mid];
	        		#else
	        		dx = tny.dlist2[mid];
	        		#endif
	        		dy = tnx.dlist2[mid];
	        	}
	        	#endif
	        	#ifdef USE_FECTH
				_mm_prefetch(dx, _MM_HINT_T0);
				_mm_prefetch(dy, _MM_HINT_T0);
				#endif
        		bGB = true;
        		break;
        	}

        
        LENGTH res = infinity;
		
        //
		if(!bGB){
			hopsize = vtree[lca].posSizes[0];
			//g_avgw += hopsize;
			//if(hopsize > g_mw)
		    //	g_mw = hopsize;
			hop = vtree[lca].positions[0];
			dx = tnx.dlist2[0];dy = tny.dlist2[0];
			#ifdef USE_FECTH
			_mm_prefetch(dx, _MM_HINT_T0);
			_mm_prefetch(dy, _MM_HINT_T0);
			_mm_prefetch(hop, _MM_HINT_T0);
			#endif
            for (int i = 0; i < hopsize; i++){
            	int j = hop[i];
            	
            	if(j >= height[x])
            		cout<< j;
            	if(j >= height[y])
            		cout<<j;
            	
                LENGTH tmp = dx[j] + dy[j];
                if (res > tmp)
                    res = tmp;
            }
		} else {

		    #ifdef USE_BIT
		    //g_avgw += hopsize;
		    //if(hopsize > g_mw)
		   // 	g_mw = hopsize;
		    int num_bytes = bv->NumByte();

		    int jj = 0; int tt = 0;
		    for(int bb = 0; bb < num_bytes; bb++){
                uint32_t word = bv->GetByte(bb);
		        while(word){
		            int ind = jj + 31 - POPCNT32(P(word));
		            #ifdef DX3
		            LENGTH tmp = dx[tt] + dy[ind];
		            ++tt;
		            #else
		            LENGTH tmp = dx[ind] + dy[ind];
		            #endif
                    word = E(word);

                    if(res > tmp)
                        res = tmp;

		        }
		        jj += 8;
		    }
		    #endif
		    #ifdef USE_BYTE
		   // g_avgw += hopsize;
		   // if(hopsize > g_mw)
		    //	g_mw = hopsize;
		    int local_pos = 0;

		    for(int i = 0; i<hopsize; i++){
		    	local_pos += bytes[i];
		    	#ifdef DX3
		    	LENGTH tmp = dx[i] + dy[local_pos];
		    	#else
		        LENGTH tmp = dx[local_pos] + dy[local_pos];
		        #endif
		        //cout<<local_pos<<'\t'<<tmp<<endl;
                if(res > tmp){
		            res = tmp;
		        }
		    }
		    #endif
		    #ifdef NO_PRUNE
			hopsize = vtree[lca].posSizes[0];
			//g_avgw += hopsize;
			//if(hopsize > g_mw)
		    //	g_mw = hopsize;
		    for(int ii = 0; ii < hopsize; ++ii)
		    {
		    	if(ii<0)
		    		cout<<ii;
		    	LENGTH tmp = dx[ii] + dy[ii];
		    	if(res>tmp)
		    		res=tmp;
		    }
		    #endif
		}
        return res;

	}
}
	int queryCnt = 1000000;
	int prune = 0;

void subq()
{
	
	int * vec1 = (int *)malloc(sizeof(int)*queryCnt);
	int * vec2 = (int *)malloc(sizeof(int)*queryCnt);

	srand((unsigned)time(NULL));
	for(int jj = 0; jj < queryCnt; ++jj)
	{
		vec1[jj] =  rand() / (double)RAND_MAX*iVNum;
		vec2[jj] =  rand() / (double)RAND_MAX*iVNum;
	}
	/*
	queryCnt = iVNum/2;
	vector<int> vTemp(iVNum);
	for(int ii = 0; ii < iVNum; ++ii)
		vTemp[ii] = ii;
	random_shuffle(vTemp.begin(), vTemp.end());
	int * vec1 = (int *)malloc(sizeof(int)*queryCnt);
	int * vec2 = (int *)malloc(sizeof(int)*queryCnt);
	for(int ii = 0, jj = 0; jj < iVNum && ii < queryCnt; ++ii, jj+=2)
	{
		vec1[ii] = vTemp[jj];
		vec2[ii] = vTemp[jj+1];
	}
	*/
	//g_mw = 0;
	timer tm;

	size_t distance = 0;

	tm.restart();

	if(prune)
	{
		for (int iii = 0; iii < queryCnt; ++iii)
		{
			#ifdef REPEAT
			int tmp = 0;
			for(int tt = 0; tt < REPEAT; ++tt)
		    	tmp += P2H(vec1[iii], vec2[iii]);
		    distance += tmp / REPEAT;
		    #else
		    distance += P2H(vec1[iii], vec2[iii]);
		    #endif
		}
	}
	else
	{
		
		for (int iii = 0; iii < queryCnt; ++iii)
		{
			#ifdef REPEAT
			int tmp = 0;
			for(int tt = 0; tt < REPEAT; ++tt)
		    	tmp += H2H(vec1[iii], vec2[iii]);
		    distance += tmp / REPEAT;
		    #else
			distance += H2H(vec1[iii], vec2[iii]);
			#endif
		}

	}
	#ifdef REPEAT
	cout << tm.getTime() / queryCnt / REPEAT<< '\t';
	#else
		cout << tm.getTime() / queryCnt << '\t';
	#endif
	cout << distance * 1.0 / queryCnt << '\t';
    cout << g_avgw * 1.0 / queryCnt << '\t' << g_mw<<endl;
    free(vec1);free(vec2);
}
int main(int argc, char* argv[])
{
	string dataset("NY");

	int method = 3;

	bool bRoLCA = true;
	queryCnt = 1000000;
	bheight_control = 10;
	prune = 0;

	int ii = 6;

	if (argc>ii)
		queryCnt = atoi(argv[ii]);

	--ii;
	if(argc>ii)
		prune = atoi(argv[ii]);

	--ii;
	if (argc>ii)
		bheight_control = atoi(argv[ii]);

	--ii;
	if (argc > ii)
		bRoLCA = atoi(argv[ii]);

	--ii;
	if (argc>ii)
		method = atoi(argv[ii]);

	--ii;
	if (argc>ii)
		dataset = argv[ii];
	string strMethod("");
	#ifdef USE_FECTH
	strMethod += "fetch\t";
	#else
	strMethod += "nofetch\t";
	#endif
	#ifdef USE_BIT
	strMethod += "bitp";
	#endif
	#ifdef USE_BYTE
	strMethod += "byte";
	#endif
	#ifdef NO_PRUNE
	strMethod += "nop";
	#endif
	#ifdef DX3
	strMethod +="3";
	#endif

	cout << strMethod<<'\t'<<dataset << '\t' ;//<< prune << "\t";//<<strMethod<<endl;

	char indexName[100];
	sprintf(indexName, "%s_%d-%d-%d.pindex", dataset.c_str(), method, bRoLCA, bheight_control);

	char queryName[100];
	sprintf(queryName, "%s.query", dataset.c_str());

	readIndex(indexName);

	LOG2 = (int*)malloc(sizeof(int) * (iVNum * 2 + 10));
	LOGD = (int*)malloc(sizeof(int) * (iVNum * 2 + 10));
	int k = 0, j = 1;
	for (int i = 0; i < iVNum * 2 + 10; i++)
	{
		if (i > j * 2)
		{
			j *= 2;
			k++;
		}
		LOG2[i] = k;
		LOGD[i] = j;
	}

/*
	const int iNumOfThread = 8;
	vector<thread> vThreads;
	for(int tt = 0; tt < iNumOfThread; ++tt)
		vThreads.push_back(thread(subq));
	for(int tt = 0; tt < iNumOfThread; ++tt)
		vThreads[tt].join();
*/
	subq();
	freeall();
	return 0;
}