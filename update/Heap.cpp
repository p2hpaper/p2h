#include "Heap.h"

#include <cstring>
//MinHeap::~MinHeap()
//{
//}

void MinHeap::initial(int _heapsize,int _indexSize)
{
	heapsize=_heapsize;
	minheap.resize(heapsize);
	if(_indexSize < 1)
		_indexSize = heapsize;
	index.resize(_indexSize);
	count=0;
	NO=_heapsize;
}


void MinHeap::Exchange(int i, int j)
{
	index[minheap[i].second] = j;
	index[minheap[j].second] = i;

	swap(minheap[i],minheap[j]);
}


void MinHeap::MinHeapIfy(int i)
{
	int l=left(i),r=right(i);
	int smallest;
	if(l<count && minheap[l]<minheap[i])
		smallest=l;
	else
		smallest=i;

	if(r<count && minheap[r]<minheap[smallest])
		smallest=r;

	if(smallest!=i)
	{
		Exchange(i,smallest);
		MinHeapIfy(smallest);
	}
}


void MinHeap::DecreaseKey(int i, int key)
{
	minheap[i].first = key;
	int p;
	while(i>0 && minheap[i]<minheap[p=parent(i)])
	{
		Exchange(i,p);
		i=p;
	}
}


void MinHeap::IncreaseKey(int i, int key)
{
	minheap[i].first=key;
	MinHeapIfy(i);
}


int MinHeap::Pop(int & vid)
{
	vid = minheap[0].second;
	index[vid]=NO;
	int key = minheap[0].first;
	count--;
	if(count==0)
		return key;
	minheap[0]=minheap[count];
	index[minheap[0].second] = 0;
	MinHeapIfy(0);

	return key;
}


void MinHeap::Push(int key, int vid)
{
	minheap[count].second=vid;
	index[vid] = count;
	count++;
	DecreaseKey(count-1,key);
}


void MinHeap::Reset()
{
	count=0;
	++NO;
}




MinHeap2::~MinHeap2()
{
	//cout<<minheap.size()<<endl;
	
	vector<pair<KEY_TYPE, int> > ().swap(minheap);
	//vector<int>().swap(index);
	
	free(index);
	index=NULL;
}

void MinHeap2::initial(int _heapsize)
{
	heapsize=_heapsize;
	minheap.resize(heapsize);

	count=0;
	NO=_heapsize;
	//index.clear();
	//index.resize(heapsize,NO);
	index=(int *)malloc(sizeof(int)*heapsize*2);
	//memset(index, NO, heapsize);
	for(int ii = 0; ii < heapsize; ++ii)
		index[ii] = NO;
}

void MinHeap2::Exchange(int i, int j)
{
	index[minheap[i].second] = j;
	index[minheap[j].second] = i;

	swap(minheap[i],minheap[j]);
}


void MinHeap2::MinHeapIfy(int i)
{
	int l=left(i),r=right(i);
	int smallest;
	if(l<count && minheap[l]<minheap[i])
		smallest=l;
	else
		smallest=i;

	if(r<count && minheap[r]<minheap[smallest])
		smallest=r;

	if(smallest!=i)
	{
		Exchange(i,smallest);
		MinHeapIfy(smallest);
	}
}


void MinHeap2::DecreaseKey(int i, KEY_TYPE key)
{
	minheap[i].first = key;
	int p;
	while(i>0 && minheap[i]<minheap[p=parent(i)])
	{
		Exchange(i,p);
		i=p;
	}
}


void MinHeap2::IncreaseKey(int i, KEY_TYPE key)
{
	minheap[i].first=key;
	MinHeapIfy(i);
}


KEY_TYPE MinHeap2::Pop(int & vid)
{
	vid = minheap[0].second;
	assert(vid < heapsize);
	index[vid]=NO;
	KEY_TYPE key = minheap[0].first;
	count--;
	if(count==0)
		return key;
	minheap[0]=minheap[count];
	assert(minheap[0].second < heapsize);
	index[minheap[0].second] = 0;
	MinHeapIfy(0);

	return key;
}


void MinHeap2::Push(KEY_TYPE key, int vid)
{
	if(index[vid] < heapsize)
		cout << index[vid] << '\t' <<vid<<endl;
	minheap[count].second=vid;
	index[vid] = count;
	count++;
	DecreaseKey(count-1,key);
	
}


void MinHeap2::Reset()
{
	count=0;
	++NO;
}


vector<int> * KEY_TYPE::s_pi=NULL;





/*
int MinHeapForVirtualDist::GetIndex(int vid, int sid)
{
	VS vs(vid, sid);
	map<VS, int>::iterator iter = mapIndex.find(vs);
	if(iter == mapIndex.end())
		return NOT_FOUND;
	else
		return iter->second;
}

LENGTH MinHeapForVirtualDist::GetLen(int index)
{
	return minheap[index].vdist;
}

*/