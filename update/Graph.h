#ifndef GRAPH_H
#define GRAPH_H

#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include "MyType.h"

using namespace std;

const int Control_DegreeBound = 60;
const bool Control_CheckCreateEdge = false;

enum NODE_TYPE{
	TYPE_Cut = 3,
	//TYPE_D2 = 2,//leaf
	TYPE_D1 = 1,//single child
	TYPE_C2 = 0,//lca
};


class TRIPLE
{
public:
	TRIPLE(int _vid, int _eid1, int _eid2):vid(_vid),eid1(_eid1),eid2(_eid2){}
	int vid;
	int eid1;
	int eid2;

	bool operator<(TRIPLE const & other) const
	{
		if(vid < other.vid)
			return true;
		else
			return false;
	}
};



//const vector<int> vZeors(64,0);
class VORDER
{
public:
	int vid;
	//int level;
	int degree;

	bool operator<(VORDER const & other) const
	{
		/*
		if (level < other.level)
		return true;
		if (level > other.level)
		return false;
		*/
		//equal level
		return (degree < other.degree);
	}
};

typedef pair<int, LENGTH> DISTANCE_ENTRY_TYPE;


struct comDENTRY
{
	bool operator() (DISTANCE_ENTRY_TYPE a, DISTANCE_ENTRY_TYPE b) { return a.first < b.first; }
};

class DISTANCE_LIST
{
public:
	vector<DISTANCE_ENTRY_TYPE > vDistance;

	int size()
	{
		return vDistance.size();
	}

	void resize(int n)
	{
		vDistance.resize(n);
	}

	void push_back(int id, LENGTH dist)
	{
		vDistance.push_back(DISTANCE_ENTRY_TYPE (id, dist));
	}

	void sortByID()
	{
		comDENTRY myComDENTRY;
		sort(vDistance.begin(), vDistance.end(), myComDENTRY);
	}

	void print(FILE * fout)
	{
		int iDLSize = vDistance.size();
		fwrite(&iDLSize, sizeof(int), 1, fout);
		//printDL(TG, vid, fout);
		for(int ii = 0; ii < iDLSize; ++ii)
		{
			fwrite(&vDistance[ii].first, sizeof(int), 1, fout);
			fwrite(&vDistance[ii].second, sizeof(LENGTH), 1, fout);
		}
	}
};

//typedef vector<LENGTH> DISTANCE_LIST_TYPE;
typedef DISTANCE_LIST DISTANCE_LIST_TYPE;

class Node{
public:
static double s_outputThredhold;
	vector<int> neighbors, pos;
	map<int/*good branch id*/, vector<int> /*pos*/> mPos;
	vector<LENGTH> vToAncestor;
	vector<LENGTH> VL;
	vector<int> children;
	int height; int branch_height;
	int parent;
	//int uniqueVertex;
	int type;
	int iCountDescendant;
	int iTimeStamp;

	Node(){
		neighbors.clear();
		neighbors.resize(0);
		//neighbors.reserve(0);
		VL.clear();
		VL.resize(0);
		//VL.reserve(0);
		pos.clear();

		children.clear();
		//the intial value are special setting, should not be set to other values for initial
		parent = -1;
		//uniqueVertex = -1;
		height = 0;
		iCountDescendant = 1;
		type = 0;
		iTimeStamp = -1;
	}

	void clear()
	{
		neighbors.clear();
		neighbors.resize(0);
		VL.clear();
		VL.resize(0);
		pos.clear();
		pos.resize(0);

		children.clear();
		children.resize(0);
		parent = -1;
		height = 0;
		iCountDescendant = 1;
		type = 0;
		iTimeStamp = -1;
	}

	bool bLeaf(const vector<int> & pi)
	{
		for (int ii = 0; ii < neighbors.size(); ++ii)
		{
			if (pi[neighbors[ii]] != -1)
				return false;
		}
		return true;
	}

	void eraseAChild(int child)
	{
		vector<int>::iterator iter = find(children.begin(), children.end(), child);
		if (iter != children.end())
			children.erase(iter);
	}

	void replaceChild(int child, int newChild)
	{
		vector<int>::iterator iter = find(children.begin(), children.end(), child);
		if (iter != children.end())
		{
			*iter = newChild;
		}
	}

	void printDL(vector<Node> & TG, int vid, FILE * fout)
	{
		if(TG[vid].vToAncestor.size() == 0)
			return;
		int parent = TG[vid].parent;
		bool bFirst = true;
		while(parent != -1)
		{
			if(bFirst)
			{
				bFirst = !(TG[parent].children.size() > 1);
				int height = TG[parent].height;
				fwrite(&TG[vid].vToAncestor[height], sizeof(LENGTH), 1, fout);
			}
			
			if(TG[parent].children.size() > 1 && !bFirst)
			{
				int height = TG[parent].height;
				fwrite(&TG[vid].vToAncestor[height], sizeof(LENGTH), 1, fout);
			}
			parent = TG[parent].parent;
		}
	}

	void printToAncestor(FILE * fout)
	{
		int iDLSize = vToAncestor.size();
		fwrite(&iDLSize, sizeof(int), 1, fout);
		//printDL(TG, vid, fout);
		fwrite(vToAncestor.data(), sizeof(LENGTH), iDLSize, fout);
	}

	void printmPosOld(FILE * fout, map<int, double> & mBidBetter)
	{
		//map<int/*branch id*/,vector<int>/*pos*/> mPos;
		int iSize = 0;
		vector<map<int, vector<int> >::iterator> vtemp;
		for(map<int, vector<int> >::iterator iter = mPos.begin(); iter != mPos.end(); ++iter)
		{
			int bid = iter->first;
			if(mBidBetter[bid] < s_outputThredhold)
			{
				++iSize;
				vtemp.push_back(iter);
			}
		}
		
		fwrite(&iSize, sizeof(int), 1, fout);
		if(iSize < 1)
			return;
		for(int ii = 0; ii < iSize; ++ii)
		{
			map<int, vector<int> >::iterator iter = vtemp[ii];
			int bid = iter->first;

			vector<int> & vpos = iter->second;
			fwrite(&bid, sizeof(int), 1, fout);
			int size = vpos.size();
			fwrite(&size, sizeof(int), 1, fout);
			if(size>0)
			{
				sort(vpos.begin(), vpos.end());
				fwrite(vpos.data(), sizeof(int), size, fout);
			}
		}
	}

	void printmPos(FILE * fout, map<int, double> & mBidBetter)
	{
		//map<int/*branch id*/,vector<int>/*pos*/> mPos;
		int iSize = mPos.size();
		fwrite(&iSize, sizeof(int), 1, fout);
		if(iSize < 1)
			return;
		for(map<int, vector<int> >::iterator iter = mPos.begin(); iter != mPos.end(); ++iter)
		{
			int bid = iter->first;

			vector<int> & vpos = iter->second;
			fwrite(&bid, sizeof(int), 1, fout);
			int size = vpos.size();
			fwrite(&size, sizeof(int), 1, fout);
			if(size>0)
			{
				sort(vpos.begin(), vpos.end());
				fwrite(vpos.data(), sizeof(int), size, fout);
			}
		}
	}

	void printMe(vector<Node> & TG, int vid, FILE * fout)
	{
		fwrite(&vid, sizeof(int), 1, fout);

		fwrite(&type, sizeof(int), 1, fout);
		bool bFake = false;

		switch (type)
		{
		case TYPE_D1:
		{
			fwrite(&height, sizeof(int), 1, fout);
			/*
			int iPosSize = pos.size();
			fwrite(&iPosSize, sizeof(int), 1, fout);
			fwrite(pos.data(), sizeof(int), iPosSize, fout);
			*/
		}
		break;
		case TYPE_C2:
		{
			//fwrite(&branch_height, sizeof(int), 1, fout);
			fwrite(&height, sizeof(int), 1, fout);
				int iPosSize = pos.size();
				fwrite(&iPosSize, sizeof(int), 1, fout);
				sort(pos.begin(), pos.end());
				fwrite(pos.data(), sizeof(int), iPosSize, fout);


		}
		break;
		default:
			assert(false);
			break;
		}

		printToAncestor(fout);

	}
};



typedef Node TreeNode;

class Vertex
{
public:
	Vertex() :degree(0.0){}
	int id;
	set<int> edge_list;

	void insertANeighborEdge(int eid)
	{
		edge_list.insert(eid);
	}

	void deleteANeighborEdge(int eid)
	{
		//assert(edge_list.find(eid) != edge_list.end());
		if(edge_list.find(eid) == edge_list.end())
			cout<<"error";
		edge_list.erase(eid);
	}

	LENGTH coordinate[2];

	double degree;


};
typedef vector<Vertex> VertexList;

class Edge{
public:
	Edge() :mark(0){}
	int id;
	LENGTH length;

	int start;
	int end;

	int mark;

	int adjVertex(int i)
	{
		if (i != start  && i != end)
		{
			cout << "vertex " << i << " not belong to edge " << id << endl;
			throw;
		}
		return i == start ? end : start;
	}
	//edge a > edge b: a.upperBound>b.upperbound;
	//bool operator > (Edge const & other) const;
	//bool operator < (Edge const & other) const;
	//bool operator >= (Edge const & other) const;
	//bool operator == (Edge const & other) const;
};

struct comVV{
	bool operator()(const pair<int, int> &v1, const pair<int, int> &v2)
	{
		if (v1.first < v2.first)
			return true;
		if (v1.first > v2.first)
			return false;
		return v1.second < v2.second;
	}
};

class EdgeList
{
public:
	vector<Edge> edges;
	map<pair<int, int>, int, comVV> mapVToE;

	int g_iCount;

	EdgeList() :g_iCount(0){}

	void clear()
	{
		edges.clear();
		edges.resize(0);

		mapVToE.clear();
	}

	void copy(const EdgeList & other)
	{
		edges = *new vector<Edge>(other.edges);
		mapVToE = *new map<pair<int, int>, int, comVV>(other.mapVToE);
	}

	void component(int vid, vector<Vertex> & vg)
	{
		set<int> snew;
		set<int> component;
		snew.insert(vid);

		while (snew.size()>0)
		{
			component.insert(snew.begin(), snew.end());
			set<int> temp;
			for (set<int>::iterator iter = snew.begin(); iter != snew.end(); ++iter)
			{
				for (set<int>::iterator eIter = vg[*iter].edge_list.begin(); eIter != vg[*iter].edge_list.end(); ++eIter)
				{
					int x = edges[*eIter].adjVertex(*iter);
					if (component.find(x) == component.end())
						temp.insert(x);
				}
			}
			snew = temp;
		}
		cout << component.size() << endl;
	}

	bool insertAnEdge(Edge & eNew)
	{
		if (eNew.start == eNew.end)
		{
			cout << "Warning: start=end when create edges" << endl;
			return false;
		}
		int eid;
		if (existEdge(eNew.start, eNew.end, eid))
		{
			if (edges[eid].length > eNew.length)
				edges[eid].length = eNew.length;
			return false;
		}
		else
		{
			pair<int, int> vv(eNew.start, eNew.end);
			mapVToE[vv] = edges.size();
			edges.push_back(eNew);
			return true;
		}
	}

	bool existEdge(int v1, int v2, int & eid)
	{
		pair<int, int> vv(v1, v2);
		map<pair<int, int>, int, comVV>::iterator iter = mapVToE.find(vv);
		if (iter != mapVToE.end())
		{
			eid = iter->second;
			assert(edges[eid].mark != -1);
			return true;
		}

		pair<int, int> vv2(v2, v1);
		iter = mapVToE.find(vv2);
		if (iter != mapVToE.end())
		{
			eid = iter->second;
			assert(edges[eid].mark != -1);
			return true;
		}

		return false;
	}

	void removeAnEdge(int eid)
	{
		edges[eid].mark = -1;
		pair<int, int> temp;
		temp.first = edges[eid].start;
		temp.second = edges[eid].end;

		map<pair<int, int>, int, comVV>::iterator iter = mapVToE.find(temp);
		if(iter == mapVToE.end())
		{
			temp.first = edges[eid].end;
			temp.second = edges[eid].start;
			iter = mapVToE.find(temp);
		}
		if(iter != mapVToE.end())
			mapVToE.erase(mapVToE.find(temp));

	}
};



class VertexInH :public Vertex
{
public:
	int level;
	int index;
	bool deleted;

	int mark;
	VertexInH(){}
	VertexInH(const Vertex & other) :Vertex(other), deleted(false), level(0), mark(0){}

	int degreeInc(vector<VertexInH> & vh, EdgeList & eh)
	{
		int degree = edge_list.size();

		if (degree > Control_DegreeBound)
			return degree;//*(degree-3)/2;
		int iResult = degree*(degree - 2);
		int iCount = 0;

		for (set<int>::iterator ii = edge_list.begin(); ii != edge_list.end(); ++ii)
		{
			int start = eh.edges[*ii].adjVertex(id);
			for (set<int>::iterator jj = ii; jj != edge_list.end(); ++jj)
			{
				if (jj == ii)
					continue;

				int end = eh.edges[*jj].adjVertex(id);
				int temp;
				if (eh.existEdge(start, end, temp))
					++iCount;
				else if (eh.existEdge(end, start, temp))
					++iCount;
			}
		}
		return iResult - 2 * iCount;

	}

	int degreeInc2(vector<VertexInH> & vh, EdgeList & eh)
	{
		int degree = edge_list.size();


		if (degree > Control_DegreeBound)
			return 2 * degree - 1;//*(degree-3)/2;

		int iMax = -1 * vh.size();
		for (set<int>::iterator ii = edge_list.begin(); ii != edge_list.end(); ++ii)
		{
			int iCount = 0;
			int start = eh.edges[*ii].adjVertex(id);
			for (set<int>::iterator jj = edge_list.begin(); jj != edge_list.end(); ++jj)
			{
				if (jj == ii)
					continue;

				int end = eh.edges[*jj].adjVertex(id);
				int temp;
				if (eh.existEdge(start, end, temp))
					++iCount;
			}
			int degreeLater = degree - 1 - iCount;
			if (degreeLater > iMax)
				iMax = degreeLater;
		}
		return iMax;

	}

	void removeInANeighbor(int eid, vector<VertexInH> & vh, EdgeList & eh)
	{
		int vid = eh.edges[eid].adjVertex(id);
		assert(vh[vid].deleted == false);

		//can be found this vertex in the neighborlist of its neighbor
		vh[vid].deleteANeighborEdge(eid);
	}

	void removeInAllNeighbor(vector<VertexInH> & vh, EdgeList & eh)
	{
		for (set<int>::iterator iter = edge_list.begin(); iter != edge_list.end(); ++iter)
		{
			removeInANeighbor(*iter, vh, eh);
		}
	}
	void commonNeighbor(vector<VertexInH> & vh, EdgeList & eh, int v1, int v2, vector<int> & vCommon, vector<pair<int, int> > & vEdge)
	{
		set<int> * pelist = &vh[v1].edge_list;
		int vsmall = v1, vlarge = v2;
		if (vh[v1].edge_list.size() > vh[v2].edge_list.size())
		{
			vsmall = v2; vlarge = v1;
			pelist = &vh[v2].edge_list;
		}

		for (set<int>::iterator iter = pelist->begin(); iter != pelist->end(); ++iter)
		{
			int eid = *iter;
			int vid = eh.edges[eid].adjVertex(vsmall);
			pair<int, int> pedge(eid, 0);
			if (eh.existEdge(vlarge, vid, pedge.second))
			{
				vCommon.push_back(vid);
				vEdge.push_back(pedge);
			}
		}
	}

	void removeMe(vector<VertexInH> & vh, EdgeList & eh, TreeNode & tn, bool bLeafNode, vector<set<TRIPLE> > & TabE, vector<LENGTH> & TabPhi)
	{
		deleted = true;
		//tn.uniqueVertex = id;
		tn.clear();
		vector<int> & neighbors = tn.neighbors;
		vector<LENGTH> & distances = tn.VL;

		for (set<int>::iterator iter = edge_list.begin(); iter != edge_list.end(); ++iter)
		{
			//get id of neighbor
			int vneighbor = eh.edges[*iter].adjVertex(id);

			if (vh[vneighbor].deleted == true)
			{
				cout << "vid=" << id << '\t' << "vneighbor=" << vneighbor << " edgeid=" << *iter << endl;
			}
			vh[vneighbor].deleteANeighborEdge(*iter); //removeInANeighbor

			neighbors.push_back(vneighbor);
			distances.push_back(eh.edges[*iter].length);


		}

		//create pairwise edges among neighbors
		
		int eid1,eid2;
		for (int ii = 0; ii < neighbors.size(); ++ii)
		{
			assert(eh.existEdge(id, neighbors[ii], eid1));
			TabPhi[eid1]=distances[ii];
			for (int jj = ii + 1; jj < neighbors.size(); ++jj)
			{
				int eid = 0;
				assert(eh.existEdge(id, neighbors[jj], eid2));
				TRIPLE triple0(id, eid1, eid2);

				if (eh.existEdge(neighbors[ii], neighbors[jj], eid))
				{
					TabE[eid].insert(triple0);
					//update distance
					LENGTH newLength = distances[ii] + distances[jj];
					if (newLength < eh.edges[eid].length)
						eh.edges[eid].length = newLength;
				}
				else
				{
					//create new edge
					Edge newEdge;
					newEdge.length = distances[ii] + distances[jj];
					//check whether to create edges or not, whether has smaller distance.
					//first, get common neighbor
					bool bCreate = true;

					newEdge.start = neighbors[ii];
					newEdge.end = neighbors[jj];
					newEdge.id = eh.edges.size();
					eh.insertAnEdge(newEdge);
					//update the neighbors
					vh[neighbors[ii]].edge_list.insert(newEdge.id);
					vh[neighbors[jj]].edge_list.insert(newEdge.id);
					TabE.resize(TabE.size()+1);
					TabPhi.resize(TabPhi.size()+1);
					assert(TabE.size()==eh.edges.size());
					TabE[newEdge.id].insert(triple0);
				}
			}
		}
	}

	void removeMeCut(vector<VertexInH> & vh, EdgeList & eh)
	{
		deleted = true;

		vector<int> neighbors;
		vector<LENGTH> distances;

		for (set<int>::iterator iter = edge_list.begin(); iter != edge_list.end(); ++iter)
		{
			//get id of neighbor
			int vneighbor = eh.edges[*iter].adjVertex(id);

			if (vh[vneighbor].deleted == true)
			{
				//cout << "vid=" << id << '\t' << "vneighbor=" << vneighbor << " edgeid=" << *iter << endl;
				continue;
			}
			vh[vneighbor].deleteANeighborEdge(*iter); //removeInANeighbor

			neighbors.push_back(vneighbor);
			distances.push_back(eh.edges[*iter].length);
		}

		//create pairwise edges among neighbors
		//no need to create edges among cut vertices since it has been consider in the singleSource process
		//no need to create edges if two vertices come from different subgraphs
		for (int ii = 0; ii < neighbors.size(); ++ii)
			for (int jj = ii + 1; jj < neighbors.size(); ++jj)
			{
				if (vh[neighbors[ii]].mark != vh[neighbors[jj]].mark)
					continue;
				int eid = 0;
				if (eh.existEdge(neighbors[ii], neighbors[jj], eid))
				{
					//update distance
					LENGTH newLength = distances[ii] + distances[jj];
					if (newLength < eh.edges[eid].length)
						eh.edges[eid].length = newLength;
				}
				else
				{
					//create new edge
					Edge newEdge;
					newEdge.length = distances[ii] + distances[jj];
					//check whether to create edges or not, whether has smaller distance.
					//first, get common neighbor
						newEdge.start = neighbors[ii];
						newEdge.end = neighbors[jj];
						newEdge.id = eh.edges.size();
						newEdge.mark = vh[neighbors[ii]].mark;
						eh.insertAnEdge(newEdge);

						//update the neighbors
						vh[neighbors[ii]].edge_list.insert(newEdge.id);
						vh[neighbors[jj]].edge_list.insert(newEdge.id);
		
				}
			}
		edge_list.clear();
	}
};

class Triangle
{
public:
	int A,B,C;//vertex
	int AB,BC,AC;//edge
};

#endif
