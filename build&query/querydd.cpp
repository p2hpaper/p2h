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
using namespace std;
#define USE_FECTH
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
int iii;
const LENGTH infinity = 999999999;

int bheight_control;
int *LOG2, *LOGD, **RMQIndex;
int *height;

int topk;
map<int, int> mGoodBranch;
int * vGoodBranch;
double * vProbBranch, * vBetterBranch;
double acceptThreshold = 0.5;
double outputThreshold;
map<int, double> mBidToBetter;

int g_mw;

class NODE
{
public:
	NODE()
	{
		pos0=pos1=pos2=pos3=pos4=pos5=NULL;
		mToAncestor=NULL;
		pos0Size=pos1Size=pos2Size=pos3Size=pos4Size=pos5Size=0;
	}

	int  pos0Size, pos1Size, pos2Size, pos3Size, pos4Size, pos5Size;
	int *pos0, *pos1, *pos2, *pos3, *pos4, *pos5;
	LENGTH * mToAncestor;
	int toRMQ;


	~NODE()
	{
		if(pos0)
			free(pos0);
		if(pos1)
			free(pos1);
		if(pos2)
			free(pos2);
		if(pos3)
			free(pos3);
		if(pos4)
			free(pos4);
		if(pos5)
			free(pos5);
		if(mToAncestor)
			free(mToAncestor);
	}

	void set(int mid, int size, FILE * fin)
	{
		switch(mid)
		{
			case 1:
			pos1Size = size;
			pos1 = (int*)malloc(sizeof(int)*size);
			fread(pos1, sizeof(int), size, fin);
			break;
			case 2:
			pos2Size = size;
			pos2 = (int*)malloc(sizeof(int)*size);
			fread(pos2, sizeof(int), size, fin);
			break;
			case 3:
			pos3Size = size;
			pos3 = (int*)malloc(sizeof(int)*size);
			fread(pos3, sizeof(int), size, fin);
			break;
			case 4:
			pos4Size = size;
			pos4 = (int*)malloc(sizeof(int)*size);
			fread(pos4, sizeof(int), size, fin);
			break;
			case 5:
			pos5Size = size;
			pos5 = (int*)malloc(sizeof(int)*size);
			fread(pos5, sizeof(int), size, fin);
			break;
			case 0:
			pos0Size = size;
			pos0 = (int*)malloc(sizeof(int)*size);
			fread(pos0, sizeof(int), size, fin);
			break;
			default:
			int * temp = (int*)malloc(sizeof(int)*size);
			fread(temp, sizeof(int), size, fin);
			free(temp);
			break;
		}
	}
};

NODE * vtree;
int avgw, iVNum;
int iRow;




LENGTH ** recordDx;
LENGTH ** recordDy;
int * recordSize;
int ** recordPos;

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

		tn.mToAncestor = (LENGTH *)malloc(sizeof(LENGTH)*iDLSize);
		fread(tn.mToAncestor, sizeof(LENGTH), iDLSize, fin);
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
	int jj = 0;
	for(int ii = 0; ii < topk; ++ii)
	{
		//cout<<vGoodBranch[ii]<<'\t'<<vProbBranch[ii]<<'\t'<<vBetterBranch[ii]<<endl;
		mBidToBetter[vGoodBranch[ii]] = vBetterBranch[ii];
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
	avgw = 0;
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

inline	LENGTH distanceQuery(int x, int y){
	if (x == y) return 0;

	NODE & tnx = vtree[x];
	NODE & tny = vtree[y];
	int lca = LCAQuery(tnx.toRMQ, tny.toRMQ);

	if(lca == x)
	{
		++avgw;
		return tny.mToAncestor[height[x]];
	}
	else if(lca == y)
	{
		++avgw;
		return tnx.mToAncestor[height[y]];
	}
	else {
		NODE & tnc = vtree[lca];
		int ps = tnc.pos0Size;
		int *p2 = tnc.pos0;
//return 0;
		avgw += ps;

		LENGTH *dx = tnx.mToAncestor, *dy = tny.mToAncestor;
		/*
				recordDx[iii] = dx;
		recordDy[iii] = dy;
		recordSize[iii] = ps;
		recordPos[iii] = p2;
		return 0;*/
#ifdef USE_FECTH
		_mm_prefetch(dx, _MM_HINT_T0);
		_mm_prefetch(dy, _MM_HINT_T0);
		_mm_prefetch(p2, _MM_HINT_T0);
#endif
		LENGTH res = infinity;
		//return 0;
//for(int jj = 0; jj < 100; ++jj)
		if(ps > g_mw)
			g_mw = ps;
		for (int i = 0; i < ps; i++){
			LENGTH tmp = dx[p2[i]] + dy[p2[i]];
			if (res > tmp)
				res = tmp;
		}
		return res;
	}
}
inline	LENGTH getDistInter(int x, int y){
	if (x == y) return 0;
	const NODE & tnx = vtree[x];
	const NODE & tny = vtree[y];
	int lca = LCAQuery(tnx.toRMQ, tny.toRMQ);

	if(lca == x)
	{
		++avgw;
		return tny.mToAncestor[height[x]];
	}
	else if(lca == y)
	{
		++avgw;
		return tnx.mToAncestor[height[y]];
	}
	else {
		int ps;
		int *p2;
		int *p1;
		int iSize2 =0;

		if(lca == vGoodBranch[0])
		{			
			ps = tnx.pos1Size;
			iSize2 = tny.pos1Size;
			
				p2 = tny.pos1;
				p1 = tnx.pos1;
			
		}
		else if(lca == vGoodBranch[1])
		{
			ps = tnx.pos2Size;

				ps = iSize2;
				p2 = tny.pos2;
				p1 = tnx.pos2;
		}
		else if(lca == vGoodBranch[2])
		{
			ps = tnx.pos3Size;
			iSize2 = tny.pos3Size;

				ps = iSize2;
				p2 = tny.pos3;

				p1 = tnx.pos3;
		}
		else if(lca == vGoodBranch[3])
		{
			ps = tnx.pos4Size;
			iSize2 = tny.pos4Size;

				ps = iSize2;
				p2 = tny.pos4;
				p1 = tnx.pos4;
		}
		else if(lca == vGoodBranch[4])
		{
			ps = tnx.pos5Size;
			iSize2 = tny.pos5Size;

				ps = iSize2;
				p2 = tny.pos5;
				p1 = tnx.pos5;
		}
		else
		{
			const NODE & tnc = vtree[lca];
			ps = tnc.pos0Size;
			p2 = tnc.pos0;

					avgw += ps;
		if(ps > g_mw)
			g_mw = ps;

		LENGTH *dx = tnx.mToAncestor, *dy = tny.mToAncestor;
		/*
		recordDx[iii] = dx;
		recordDy[iii] = dy;
		recordSize[iii] = ps;
		recordPos[iii] = p2;
		return 0;*/
#ifdef USE_FECTH
		_mm_prefetch(dx, _MM_HINT_T0);
		_mm_prefetch(dy, _MM_HINT_T0);
		_mm_prefetch(p2, _MM_HINT_T0);
#endif
		//return 0;
		LENGTH res = infinity;
//		for(int jj = 0; jj < 100; ++jj)
		for (int i = 0; i < ps; i++){
			LENGTH tmp = dx[p2[i]] + dy[p2[i]];
			if (res > tmp)
				res = tmp;
		}
		return res;
		}
//return 0;


		LENGTH *dx = tnx.mToAncestor, *dy = tny.mToAncestor;
		/*
		recordDx[iii] = dx;
		recordDy[iii] = dy;
		recordSize[iii] = ps;
		recordPos[iii] = p2;
		return 0;*/
#ifdef USE_FECTH
		_mm_prefetch(dx, _MM_HINT_T0);
		_mm_prefetch(dy, _MM_HINT_T0);
		_mm_prefetch(p2, _MM_HINT_T0);
#endif
		//return 0;
		LENGTH res = infinity;
		int temp = 0;
//		for(int jj = 0; jj < 100; ++jj)
		for (int i = 0, j=0; i < ps && j<iSize2; ){
			if(p1[i] < p2[j])
			{
				++i;
			}
			else if(p1[i] > p2[j])
			{
				++j;
			}
			else{
			LENGTH tmp = dx[p2[j]] + dy[p2[j]];
			++j; ++i; ++temp;
			if (res > tmp)
				res = tmp;
			}
		}
				avgw += temp;
		if(temp > g_mw)
			g_mw = temp;
		return res;
	}
}

inline	LENGTH getDist(int x, int y){
	if (x == y) return 0;
	const NODE & tnx = vtree[x];
	const NODE & tny = vtree[y];
	int lca = LCAQuery(tnx.toRMQ, tny.toRMQ);

	if(lca == x)
	{
		++avgw;
		return tny.mToAncestor[height[x]];
	}
	else if(lca == y)
	{
		++avgw;
		return tnx.mToAncestor[height[y]];
	}
	else {
		int ps;
		int *p2;

		if(lca == vGoodBranch[0])
		{			
			ps = tnx.pos1Size;
			int iSize2 = tny.pos1Size;
			
			if(ps > iSize2)
			{
				ps = iSize2;
				p2 = tny.pos1;
			}
			else
				p2 = tnx.pos1;
			
		}
		else if(lca == vGoodBranch[1])
		{
			ps = tnx.pos2Size;
			int iSize2 = tny.pos2Size;
			if(ps > iSize2)
			{
				ps = iSize2;
				p2 = tny.pos2;
			}
			else
				p2 = tnx.pos2;
		}
		else if(lca == vGoodBranch[2])
		{
			ps = tnx.pos3Size;
			int iSize2 = tny.pos3Size;
			if(ps > iSize2)
			{
				ps = iSize2;
				p2 = tny.pos3;
			}
			else
				p2 = tnx.pos3;
		}
		else if(lca == vGoodBranch[3])
		{
			ps = tnx.pos4Size;
			int iSize2 = tny.pos4Size;
			if(ps > iSize2)
			{
				ps = iSize2;
				p2 = tny.pos4;
			}
			else
				p2 = tnx.pos4;
		}
		else if(lca == vGoodBranch[4])
		{
			ps = tnx.pos5Size;
			int iSize2 = tny.pos5Size;
			if(ps > iSize2)
			{
				ps = iSize2;
				p2 = tny.pos5;
			}
			else
				p2 = tnx.pos5;
		}
		else
		{
			const NODE & tnc = vtree[lca];
			ps = tnc.pos0Size;
			p2 = tnc.pos0;
		}
//return 0;
		avgw += ps;
		if(ps > g_mw)
			g_mw = ps;

		LENGTH *dx = tnx.mToAncestor, *dy = tny.mToAncestor;
		/*
		recordDx[iii] = dx;
		recordDy[iii] = dy;
		recordSize[iii] = ps;
		recordPos[iii] = p2;
		return 0;*/
#ifdef USE_FECTH
		_mm_prefetch(dx, _MM_HINT_T0);
		_mm_prefetch(dy, _MM_HINT_T0);
		_mm_prefetch(p2, _MM_HINT_T0);
#endif
		//return 0;
		LENGTH res = infinity;
//		for(int jj = 0; jj < 100; ++jj)
		for (int i = 0; i < ps; i++){
			LENGTH tmp = dx[p2[i]] + dy[p2[i]];
			if (res > tmp)
				res = tmp;
		}
		return res;
	}
}

void prepareBatch(int iVNum, int iBatch, int iSize, string filename)
{
	int ** m1 = (int**)malloc(sizeof(int*)*iBatch);
	int ** m2 = (int**)malloc(sizeof(int*)*iBatch);

	for(int ii = 0; ii < iBatch; ++ii)
	{
		m1[ii] = (int*)malloc(sizeof(int)*iSize);
		m2[ii] = (int*)malloc(sizeof(int)*iSize);
	}

	int * vcount = (int*)malloc(sizeof(int)*iBatch);
	memset(vcount, 0, sizeof(int)*iBatch);

	const int queryCnt = 1000000;
	LENGTH iMax = 0;

	LENGTH * dis = (LENGTH *)malloc(sizeof(LENGTH)*queryCnt);
	double lall = 0;
	for(int jj = 0; jj < queryCnt; ++jj)
	{
		int v1 = rand() / (double)RAND_MAX*iVNum;
		int v2 = rand() / (double)RAND_MAX*iVNum;

		dis[jj]=distanceQuery(v1, v2);
		//if(jj < 0)
			//cout<<vec1[jj]<<'\t'<<vec2[jj]<<'\t'<<dis[jj]<<'\t';
		lall += dis[jj];
		if(dis[jj] > iMax)
			iMax = dis[jj];
	}

	cout<<"iMax="<<iMax<<"\tavgdis="<<lall/queryCnt<<endl;

	const LENGTH iMin = 1000;
	double ratio = pow(iMax*1.0/iMin, 1.0/iBatch);
	cout<<"ratio="<<ratio<<endl;
	double * bound = (double*)malloc(sizeof(double)*iBatch);
	double temp = (double)iMin;
	for(int ii = 0; ii < iBatch; ++ii)
	{
		temp *= ratio;
		bound[ii] = temp;
		cout<<bound[ii]<<'\t';
	}
	cout<<endl;

	int iNow = 0;
	const int iAll = iBatch * iSize;
	int loop = 0;
	while(iNow < iAll && loop<100000000)
	{
		++loop;
		int v1 = rand() / (double)RAND_MAX*iVNum;
		int v2 = rand() / (double)RAND_MAX*iVNum;
		LENGTH dd = distanceQuery(v1, v2);
		for(int ii = 0; ii < iBatch; ++ii)
		{
			if(dd < bound[ii])
			{
				if(vcount[ii] >= iSize)
					break;
				else{
					m1[ii][vcount[ii]] = v1;
					m2[ii][vcount[ii]] = v2;
					++vcount[ii];
					++iNow;
					break;
				}
			}
		}
	}

	//cout<<"finish preparation"<<endl;
	timer tm;
	vector<double> tPrune(iBatch);
	vector<double> tT2(iBatch);
	for(int ii = 0; ii < iBatch; ++ii)
	{
		//cout<<vcount[ii]<<'\t';
		if(vcount[ii] == 0)
		{
			cout<<"batch"<<ii<<'\t';
			continue;
		}

		int * vec1 = m1[ii];
		int * vec2 = m2[ii];

			avgw = 0;
			tm.restart();
			for(int jj = 0; jj < vcount[ii]; ++jj)
			{
				distanceQuery(vec1[jj], vec2[jj]);
			}
			tT2[ii] = tm.getTime() / vcount[ii];

			avgw = 0;
			tm.restart();
			for(int jj = 0; jj < vcount[ii]; ++jj)
			{
				getDist(vec1[jj], vec2[jj]);
			}
			tPrune[ii] = tm.getTime() / vcount[ii];
		

	}

	for(int ii = 0; ii < iBatch; ++ii)
		cout<<tPrune[ii]<<'\t'<<endl;
	cout<<endl;
	for(int ii = 0; ii < iBatch; ++ii)
		cout<<tT2[ii]<<'\t'<<endl;
	cout<<endl;
	FILE * fout = fopen(filename.c_str(), "wb");
	fwrite(&iBatch, sizeof(int), 1, fout);
	fwrite(&iSize, sizeof(int), 1, fout);
	fwrite(vcount,sizeof(int), iBatch, fout);
	fwrite(bound, sizeof(double), iBatch, fout);
	for(int ii = 0; ii < iBatch; ++ii)
	{
		fwrite(m1[ii], sizeof(int), vcount[ii], fout);
		free(m1[ii]);
		fwrite(m2[ii], sizeof(int), vcount[ii], fout);
		free(m2[ii]);
	}
	fclose(fout);

	free(m1);
	free(m2);
	free(vcount);
	free(bound);
}

void readQuery(string queryName, bool bPrune)
{
	FILE * fin = fopen(queryName.c_str(), "rb");

	int iBatch, iSize;
	fread(&iBatch, sizeof(int), 1, fin);
	fread(&iSize, sizeof(int), 1, fin);

	int * vcount = (int *)malloc(sizeof(int)*iBatch);
	fread(vcount, sizeof(int), iBatch, fin);
	double * bound = (double *)malloc(sizeof(double)*iBatch);
	fread(bound, sizeof(double), iBatch, fin);


	int ** m1 = (int**)malloc(sizeof(int*)*iBatch);
	int ** m2 = (int**)malloc(sizeof(int*)*iBatch);

	for(int ii = 0; ii < iBatch; ++ii)
	{
		m1[ii] = (int*)malloc(sizeof(int)*vcount[ii]);
		fread(m1[ii], sizeof(int), vcount[ii], fin);
		m2[ii] = (int*)malloc(sizeof(int)*vcount[ii]);
		fread(m2[ii], sizeof(int), vcount[ii], fin);
	}
	fclose(fin);

	timer tm;

	vector<double> tPrune(iBatch);
	vector<double> tT2(iBatch);
	cout<<queryName<<endl;
	if(bPrune)
		cout << "P2H\tavgw\n";
	else
		cout << "T2\tavgw\n";
	double dSum = 0;
	for(int ii = iBatch-1; ii > -1; --ii)
	{
		//cout<<vcount[ii]<<'\t';
		if(vcount[ii] == 0)
		{
			cout<<"batch"<<ii<<'\t';
			continue;
		}

		int * vec1 = m1[ii];
		int * vec2 = m2[ii];

if(!bPrune){
			avgw = 0;
			tm.restart();
			for(int jj = 0; jj < vcount[ii]; ++jj)
			{
				dSum+=distanceQuery(vec1[jj], vec2[jj]);
			}
			tT2[ii] = tm.getTime() / vcount[ii];
			cout<<1.0*avgw/vcount[ii]<<'\t';
}
else{
			avgw = 0;
			tm.restart();
			for(int jj = 0; jj < vcount[ii]; ++jj)
			{
				dSum+=getDist(vec1[jj], vec2[jj]);
			}
			tPrune[ii] = tm.getTime() / vcount[ii];
			cout<<1.0*avgw/vcount[ii]<<'\t';
		}
	}
	cout<<endl<<queryName<<'\t'<<dSum<<endl;
	if(bPrune){
	cout<<"P2H\tqtime\n";
	for(int ii = 0; ii < iBatch; ++ii)
		cout<<tPrune[ii]<<'\t';
}
else{
	cout<<"T2\tqtime\n";
	for(int ii = 0; ii < iBatch; ++ii)
		cout<<tT2[ii]<<'\t';
}
	cout<<endl;

	for(int ii = 0; ii < iBatch; ++ii)
	{
		free(m1[ii]);
		free(m2[ii]);
	}

	free(m1);
	free(m2);
	free(vcount);
	free(bound);
}

int main(int argc, char* argv[])
{
	string dataset("NY");

	int method = 3;//

	bool bRoLCA = true;
	int queryCnt = 1000000;
	bheight_control = 10;
	int prune = 0;

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

	//cout << dataset << '\t' ;//<< method << '\t';


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


readQuery(queryName, prune);
return 0;


/*
prepareBatch(iVNum, 10, 10000, queryName);
freeall();
return 0;
*/

	int * vec1 = (int *)malloc(sizeof(int)*queryCnt);
	int * vec2 = (int *)malloc(sizeof(int)*queryCnt);
	for(int jj = 0; jj < queryCnt; ++jj)
	{
		vec1[jj] =  rand() / (double)RAND_MAX*iVNum;
		vec2[jj] =  rand() / (double)RAND_MAX*iVNum;
	}
	g_mw = 0;
	timer tm;

	tm.restart();

	if(prune)
	{
		for (iii = 0; iii < queryCnt; ++iii)
		{
			//for(int jj=0;jj<2;++jj)
				getDist(vec1[iii], vec2[iii]);
		}
	}
	else
	{
		for (iii = 0; iii < queryCnt; ++iii)
		{
			//for(int jj=0;jj<2;++jj)
			distanceQuery(vec1[iii], vec2[iii]);
		}
	}
	cout << tm.getTime() / queryCnt << '\t';
	cout << avgw*1.0 / queryCnt << '\t'<<g_mw<<endl;
	
	//cout<< test.dTimeSum<<endl;
	//cout<<iCount<<'\t'<<iloop<<'\t'<<icomp<<endl;
	freeall();
	return 0;
}