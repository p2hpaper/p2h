#ifndef HEAP_H
#define HEAP_H
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include "MyType.h"
#include "Graph.h"
using namespace std;

class VS
{
public:
	VS(int _vid = 0, int _sid = 0):vid(_vid),sid(_sid){}
	int vid,sid;

	bool operator < (const VS & other) const
	{
		if(vid < other.vid)
			return true;
		else if( vid > other.vid)
			return false;
		else
			return sid < other.sid;
	}
};




class MinHeap
{
public:
	vector<pair<int, int> > minheap;//degree and vid
	vector<int> index;
	int heapsize;
	int count;
	//出栈标识符
	int NO;

	MinHeap() :heapsize(0), count(0)
	{
	}
	~MinHeap()
	{
		minheap.resize(0);
		minheap.clear();
		index.resize(0);
		index.clear();
	}
	
	inline int parent(int i)
	{
		return ((i-1)>>1);
	}
	inline int left(int i)
	{
		return (i<<1)+1;
	}
	inline int right(int i)
	{
		return (i<<1)+2;
	}
	inline bool isEmpty()
	{
		return count==0;
	}
	inline int top()
	{
		return minheap[0].first;
	}

	void initial(int _heapsize, int _indexSize = -1);

	void Exchange(int i,int j);

	void MinHeapIfy(int i);
	int Pop(int & value);
	void Push(int key, int value);
	void IncreaseKey(int i, int keyValue);
	void DecreaseKey(int i, int keyValue);
	void UpdateKey(int i, int keyValue)
	{
		int oldKey = minheap[i].first;
		if (oldKey > keyValue)
			DecreaseKey(i, keyValue);
		else if (oldKey < keyValue)
			IncreaseKey(i, keyValue);
	}
	void Reset();

	//bool Find(LENGTH keyValue,int &index,int from);
};


typedef int FIRST_ORDER;
typedef int SECOND_ORDER;
typedef int THIRD_ORDER;
class KEY_TYPE
{
public:
	FIRST_ORDER key1;
	SECOND_ORDER key2;
	THIRD_ORDER key3;


	KEY_TYPE(FIRST_ORDER _key1=0, SECOND_ORDER _key2=0, THIRD_ORDER _key3=0):key1(_key1),key2(_key2),key3(_key3){}

	bool operator < (const KEY_TYPE & other) const
	{


		if(key1 < other.key1)
			return true;
		if(other.key1 < key1)
			return false;

		if(key2 < other.key2)
			return true;
		if(other.key2 < key2)
			return false;

		if(key3 < other.key3)
			return true;
		if(other.key3 < key3)
			return false;
		return false;
	}

	bool operator > (const KEY_TYPE & other) const
	{


		if(key1 > other.key1)
			return true;
		if(other.key1 > key1)
			return false;

		if(key2 > other.key2)
			return true;
		if(other.key2 > key2)
			return false;

		if(key3 > other.key3)
			return true;
		if(other.key3 > key3)
			return false;
		return false;
	}

};

class MinHeap2
{
public:
	vector<pair<KEY_TYPE, int> > minheap;//degree and vid
	vector<int> index;
	int heapsize;
	int count;
	//出栈标识符
	int NO;

	MinHeap2() :heapsize(0), count(0)
	{
	}
	~MinHeap2();
	
	inline int parent(int i)
	{
		return ((i-1)>>1);
	}
	inline int left(int i)
	{
		return (i<<1)+1;
	}
	inline int right(int i)
	{
		return (i<<1)+2;
	}
	inline bool isEmpty()
	{
		return count==0;
	}

	void initial(int _heapsize);

	void Exchange(int i,int j);

	void MinHeapIfy(int i);
	KEY_TYPE Pop(int & value);
	void Push(KEY_TYPE key, int value);
	void IncreaseKey(int i, KEY_TYPE keyValue);
	void DecreaseKey(int i, KEY_TYPE keyValue);
	void UpdateKey(int i, KEY_TYPE keyValue)
	{
		KEY_TYPE oldKey = minheap[i].first;
		if (keyValue < oldKey)
			DecreaseKey(i, keyValue);
		else if (oldKey < keyValue)
			IncreaseKey(i, keyValue);
	}
	void Reset();

	//bool Find(LENGTH keyValue,int &index,int from);
};

enum 
{
	OUT_HEAP = -1,
	NOT_FOUND = -2,
};
#endif
