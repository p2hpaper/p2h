// h2h.cpp : Defines the entry point for the console application.
//
#include <algorithm>
#include "Readfile.h"
#include "Heap.h"
#include "lca.h"
#include "MyType.h"
#include <list>
#include <fstream>
#include <numeric>
#include <memory.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include <sys/time.h>

//#define SHOW_DETAIL
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

double file_size(const char* filename)  
{  
    struct stat statbuf;  
    stat(filename,&statbuf);  
    long size=statbuf.st_size;  
  
    return size/1000.0/1000.0;  
}

enum METHOD{
	DegreeOnly = 0,
	LB = 1,
	PotentialHeight = 2,
	NoIncDeg = 3,
	CutVersion = 4,
	COOR = 5,
};

class Tree_Decomposition{
public:
	int g_Samples;
	int g_bRoLCA;
	Tree_Decomposition(ReadFile & rf) :original_pvl(&rf.vertexList), original_pel(&rf.edgeList), pvl(NULL), pel(NULL), g_Samples(1000000),g_bRoLCA(false)
	{
	}

	void run(int method)
	{
		initial();

		switch (method)
		{
		case DegreeOnly:
			treeDecOnlyDegree();
			break;
		case PotentialHeight:
			treeDecPotentialHeight();
			break;
		case NoIncDeg:
			treeDecNoIncDeg();
			break;
		case CutVersion:
			runOneCut2();
			break;
		case COOR:
		{
		vector<int> vPriority;
		setOrder(vPriority);
			treeDecDegOrd(vPriority);
		}
			break;
		default:
			treeDecOnlyDegree();
			break;
		}

		maketree();
	}



	map<char, double> mResult;

	VertexList * pvl;
	EdgeList * pel;

	VertexList * original_pvl;
	EdgeList * original_pel;

	vector<VertexInH> vh;
	int iVNum;
	EdgeList eh;
	vector<TreeNode> TG;
	vector<int> pi;
	int order;
	int root;
	int maxHeight;

	LCA lca;

	vector<int> vGoodBranch;
	vector<double> vBetterBranch;
	vector<double> vProbBranch;

	void initial()
	{
		pvl = original_pvl;
		pel = original_pel;

		iVNum = pvl->size();
		vh.clear();
		vh.resize(0);
		//copy vl to vh
		for (int ii = 0; ii < iVNum; ++ii)
		{
			VertexInH temp(pvl->at(ii));
			vh.push_back(temp);

			//TG[ii].iReal = ii;
		}

		eh.clear();

		eh.copy(*pel);
		//dealTriangle();//iVnum is updated.

		TG.resize(iVNum);
		pi.resize(iVNum, -1);
		order = 0;
		root = 0;

		lca.initial(TG, iVNum);
	}


	void treeDecOnlyDegree()
	{
		//ofstream fout("debug.txt");
		MinHeap minheap;
		minheap.initial(iVNum);

		for (int ii = 0; ii < iVNum; ++ii)
		{
			if (!vh[ii].deleted)
				minheap.Push(vh[ii].edge_list.size(), ii);
		}

		while (minheap.count)
		{

			int vid;
			minheap.Pop(vid);

			//delete the vertex in h
			pi[vid] = order;
			++order;
			TreeNode & tn = TG[vid];

			vh[vid].removeMe(vh, eh, tn, false);

			//update the degree of its neighbors
			for (int jj = 0; jj < tn.neighbors.size(); ++jj)
			{
				int vneighbor = tn.neighbors[jj];
				int degree = vh[vneighbor].edge_list.size();

				int pos = minheap.index[vneighbor];
				minheap.UpdateKey(pos, degree);
			}
		}
	}

	void treeDecPotentialHeight()
	{
		MinHeap minheap;
		minheap.initial(iVNum);

		for (int ii = 0; ii < iVNum; ++ii)
		{
			if (!vh[ii].deleted)
				minheap.Push(vh[ii].edge_list.size(), ii);
		}


		vector<int> levels(iVNum, 1);
		while (minheap.count)
		{
			int vid;
			int currentPotentialHeight = minheap.Pop(vid);
			TreeNode & tn = TG[vid];

			//check whether need to update key
			int maxPotentialHeight = levels[vid] + tn.neighbors.size();
			for (int jj = 0; jj < tn.neighbors.size(); ++jj)
			{
				int vneighbor = tn.neighbors[jj];
				int pos = minheap.index[vneighbor];
				if (minheap.minheap[pos].first > maxPotentialHeight)
					maxPotentialHeight = minheap.minheap[pos].first;
			}
			if (maxPotentialHeight > currentPotentialHeight)
			{
				minheap.Push(maxPotentialHeight, vid);
				continue;
			}

			//delete the vertex in h
			pi[vid] = order;
			++order;
			vh[vid].removeMe(vh, eh, tn, levels[vid] == 1);

			//update the potential height of its neighbors
			for (int jj = 0; jj < tn.neighbors.size(); ++jj)
			{
				int vneighbor = tn.neighbors[jj];
				int degree = vh[vneighbor].edge_list.size();

				levels[vneighbor] = Max(levels[vid] + 1, levels[vneighbor]);


				int pos = minheap.index[vneighbor];
				minheap.UpdateKey(pos, degree + levels[vneighbor]);
			}
		}
	}

	void treeDecNoIncDeg()
	{
		MinHeap2 minheap;
		minheap.initial(iVNum);

		for (int ii = 0; ii < iVNum; ++ii)
		{
			if (!vh[ii].deleted)
			{
				KEY_TYPE key(vh[ii].edge_list.size(), /*vh[ii].degreeInc(vh, eh),*/ vh[ii].edge_list.size());
				minheap.Push(key, ii);
			}
		}

		vector<int> levels(iVNum, 1);
		while (minheap.count)
		{
			int vid;
			KEY_TYPE currentKey = minheap.Pop(vid);
			TreeNode & tn = TG[vid];

			//check whether need to update key
			//KEY_TYPE maxPotentialHeight(currentKey.key1, vh[vid].degreeInc(vh,eh));

			//if(currentKey < maxPotentialHeight)
			//{
			//	minheap.Push(maxPotentialHeight, vid);
			//	continue;
			//}

			//delete the vertex in h
			pi[vid] = order;
			++order;
			vh[vid].removeMe(vh, eh, tn, levels[vid] == 1);


			//update the potential height of its neighbors
			for (int jj = 0; jj < tn.neighbors.size(); ++jj)
			{
				int vneighbor = tn.neighbors[jj];
				int degree = vh[vneighbor].edge_list.size();

				levels[vneighbor] = Max(levels[vid] + 1, levels[vneighbor]);


				int pos = minheap.index[vneighbor];
				KEY_TYPE temp(degree + levels[vneighbor],  /*vh[vneighbor].degreeInc(vh, eh),*/ degree);
				minheap.UpdateKey(pos, temp);
			}
		}
	}

	void setOrder(vector<int> & vPriority)
	{
		vPriority.resize(iVNum,-1);
		//normalize the coordinate.
		LENGTH minx = pvl->at(0).coordinate[0], miny = pvl->at(0).coordinate[1];
		LENGTH maxx = minx, maxy = miny;
		for(int ii = 1; ii < iVNum; ++ii)
		{
			LENGTH x = pvl->at(ii).coordinate[0], y = pvl->at(ii).coordinate[1];
			if(x < minx)
				minx = x;
			else if(x > maxx)
				maxx = x;

			if(y < miny)
				miny = y;
			else if(y > maxy)
				maxy = y;
		}
		double cx = (maxx + minx) / 2.0, cy = (maxy + miny) / 2.0;
		double scale = (maxy - miny) / (maxx - minx);
		scale *= scale;
		for(int ii = 0; ii < iVNum; ++ii)
		{
			LENGTH x = pvl->at(ii).coordinate[0], y = pvl->at(ii).coordinate[1];
			vPriority[ii] = -(int)( (x - cx)*(x - cx)*scale + (y - cy)*(y - cy) );
		}

	}

	void treeDecDegOrd(vector<int> &vPriority)
	{
		
		MinHeap2 minheap;
		minheap.initial(iVNum);

		for (int ii = 0; ii < iVNum; ++ii)
		{
			if (!vh[ii].deleted)
			{
				int iPriotiry = vPriority[ii];
				if(iPriotiry == -1)
				{
					KEY_TYPE key(vh[ii].edge_list.size(), vh[ii].edge_list.size());
					minheap.Push(key, ii);
				}
				else
				{
					KEY_TYPE key(vh[ii].edge_list.size(), iPriotiry);
					minheap.Push(key, ii);
				}
			}
		}


		while (minheap.count)
		{
			int vid;
			KEY_TYPE currentKey = minheap.Pop(vid);
			TreeNode & tn = TG[vid];

			//if(minheap.count < 72)
				//cout<<vid<<","<<vh[vid].edge_list.size()<<endl;

			//check whether need to update key
			//KEY_TYPE maxPotentialHeight(currentKey.key1, vh[vid].degreeInc(vh,eh));

			//if(currentKey < maxPotentialHeight)
			//{
			//	minheap.Push(maxPotentialHeight, vid);
			//	continue;
			//}

			//delete the vertex in h
			pi[vid] = order;
			++order;
			vh[vid].removeMe(vh, eh, tn, false);


			//update the potential height of its neighbors
			for (int jj = 0; jj < tn.neighbors.size(); ++jj)
			{
				int vneighbor = tn.neighbors[jj];
				int degree = vh[vneighbor].edge_list.size();

				int pos = minheap.index[vneighbor];
				KEY_TYPE temp(degree,vPriority[vneighbor]);
				//if(vPriority[vneighbor] == -1)
					minheap.UpdateKey(pos, temp);
			}
		}
	}

	void treeDecOrdDeg(vector<int> & vPriority)
	{
		MinHeap2 minheap;
		minheap.initial(iVNum);

		for (int ii = 0; ii < iVNum; ++ii)
		{
			if (!vh[ii].deleted)
			{
				int iPriotiry = vPriority[ii];
				if(iPriotiry == -1)
				{
					KEY_TYPE key(vh[ii].edge_list.size(), /*vh[ii].degreeInc(vh, eh),*/ vh[ii].edge_list.size());
					minheap.Push(key, ii);
				}
				else
				{
					KEY_TYPE key(iPriotiry, /*vh[ii].degreeInc(vh, eh),*/ vh[ii].edge_list.size());
					minheap.Push(key, ii);
				}
			}
		}

		vector<int> levels(iVNum, 1);
		while (minheap.count)
		{
			int vid;
			KEY_TYPE currentKey = minheap.Pop(vid);
			TreeNode & tn = TG[vid];

			//if(minheap.count < 72)
				//cout<<vid<<","<<vh[vid].edge_list.size()<<endl;

			//check whether need to update key
			//KEY_TYPE maxPotentialHeight(currentKey.key1, vh[vid].degreeInc(vh,eh));

			//if(currentKey < maxPotentialHeight)
			//{
			//	minheap.Push(maxPotentialHeight, vid);
			//	continue;
			//}

			//delete the vertex in h
			pi[vid] = order;
			++order;
			vh[vid].removeMe(vh, eh, tn, levels[vid] == 1);


			//update the potential height of its neighbors
			for (int jj = 0; jj < tn.neighbors.size(); ++jj)
			{
				int vneighbor = tn.neighbors[jj];
				int degree = vh[vneighbor].edge_list.size();

				levels[vneighbor] = Max(levels[vid] + 1, levels[vneighbor]);


				int pos = minheap.index[vneighbor];
				KEY_TYPE temp(degree + levels[vneighbor],  /*vh[vneighbor].degreeInc(vh, eh),*/ degree);
				if(vPriority[vneighbor] == -1)
					minheap.UpdateKey(pos, temp);
			}
		}
	}

	void maketree()
	{
		int width = 0;
		//start from the last deleted vertex, which is the root
		vector<int> ord(iVNum, 0);
		for (int ii = 0; ii < iVNum; ++ii)
		{
			ord[pi[ii]] = ii; //ord[jj] = ii : the vertex ii are the jj=pi[ii]-th deleted vertex
		}
		for (int jj = iVNum - 1; jj >= 0; --jj)
		{
			int v = ord[jj];

			TreeNode & tn = TG[v];

			int iNumNeighbor = tn.neighbors.size();
			if (iNumNeighbor > width)
				width = iNumNeighbor;

			if (iNumNeighbor > 0)
			{
				int least = iVNum;
				int u = 0;
				for (int jj = 0; jj < iNumNeighbor; ++jj)
				{
					int vneighbor = tn.neighbors[jj];
					if (pi[vneighbor] < least)
					{
						least = pi[vneighbor];
						u = vneighbor;
					}
				}
				tn.parent = u;
				TG[u].children.push_back(v);
				TG[v].height = TG[u].height + 1;
			}
			else
			{
				root = v;
				TG[root].height = 0;
#ifdef SHOW_DETAIL
				cout << "root = " << root << endl;
#else
				cout << root << '\t';
#endif
			}
		}

		int leaves = 0;
		int branches=0;
		maxHeight = 0;
		for (int v = 0; v < iVNum; ++v)
			if (TG[v].children.size() == 0)
			{
				++leaves;
				if (TG[v].height > maxHeight)
					maxHeight = TG[v].height;
				TG[v].type=TYPE_D1;
			}
			else if(TG[v].children.size() == 1)
			{
				TG[v].type=TYPE_D1;
			}
			else
			{
				++branches;
			}
#ifdef SHOW_DETAIL
		cout << "leaves = " << leaves << endl;
		cout << "height = " << maxHeight << endl;
		cout << "width = " << width << endl;
#else
		cout << "branches = " << branches << endl;
		cout << leaves << '\t' << maxHeight << '\t' << width << '\t';
#endif
	}

	void calDisAndPos(int iRoot)
	{
		//should not be called before the height is set correctly
		//calculate the distance to each ancestors

		TreeNode & tn = TG[iRoot];

		tn.vToAncestor.resize(tn.height + 1);
		tn.vToAncestor[tn.height] = 0;
		vector<int> & vNeighbor = tn.neighbors;
		vector<LENGTH> & vToNeighbor = tn.VL;
		int ancestor = tn.parent;
		while (ancestor != -1)
		{
			//calculate the distance to the ancestor
			//enumerate all pathes through the neighbors
			LENGTH minDis = (int)INFINITY;
			int posOfAncestor = TG[ancestor].height;

			for (int ii = 0; ii < vNeighbor.size(); ++ii)
			{
				int nvid = vNeighbor[ii];
				LENGTH temp = vToNeighbor[ii];
				if (posOfAncestor < TG[nvid].vToAncestor.size())
					temp += TG[nvid].vToAncestor[posOfAncestor];
				else
				{
					int posOfNvid = TG[nvid].height;
					temp += TG[ancestor].vToAncestor[posOfNvid];
				}
				if (temp < minDis)
					minDis = temp;
			}

			tn.vToAncestor[posOfAncestor] = minDis;

			ancestor = TG[ancestor].parent;
		}

		//set the pos be the height of the neighbors
		vector<int> & vPos = tn.pos;
		vPos.resize(vNeighbor.size()+1);
		for (int ii = 0; ii < vNeighbor.size(); ++ii)
			vPos[ii] = TG[vNeighbor[ii]].height;
		vPos[vNeighbor.size()] = tn.height;
		//recursive on children
		for (int ii = 0; ii < tn.children.size(); ++ii)
			calDisAndPos(tn.children[ii]);
	}


	void RoLCA(int iRoot)
	{
		TreeNode & tn = TG[iRoot];
		vector<int> & vChildren = tn.children;

		int nChildren = vChildren.size();
		if (nChildren == 0)
			return;
		if (nChildren == 1)
			RoLCA(vChildren[0]);
		else
		{			
			if(nChildren == 2)
			{
				int c1 = vChildren[0], c2 = vChildren[1];
				int c = c1;
				if(TG[c1].pos.size() > TG[c2].pos.size())
					c = c2;
				if(TG[c].pos.size() > 0)
				{
					tn.pos.resize(0);
					tn.pos.clear();
					tn.pos = TG[c].pos;

					tn.pos.resize(TG[c].pos.size()-1);
				}
			}
			else
			{
				int iMaxIndex = 0;
				int iMaxSize = 0;
				for (int ii = 0; ii < nChildren; ++ii)
				{
					int nvid = vChildren[ii];
					const int iSize = TG[nvid].neighbors.size();
					if(iSize > iMaxSize)
					{
						iMaxIndex = ii;
						iMaxSize = iSize;
					}
				}
				set<int> sTemp;
				for (int ii = 0; ii < nChildren; ++ii)
				{
					if(ii == iMaxIndex)
						continue;
					int nvid = vChildren[ii];
					const vector<int> & vn = TG[nvid].neighbors;
					const int iSize = vn.size();

					for(int jj=0; jj<iSize; ++jj)
						sTemp.insert(vn[jj]);
				}
				const int iSizeOld = tn.pos.size();
				if(tn.pos.size() > sTemp.size())
				{
					tn.pos.resize(0);
					tn.pos.clear();
					for(set<int>::iterator iter = sTemp.begin(); iter != sTemp.end(); ++iter)
						tn.pos.push_back(TG[*iter].height);
					//	tn.pos.push_back(*iter);
				}
				if(tn.pos.size() > iSizeOld)
				{
					cout << tn.pos.size()<<","<<iSizeOld<<","<<sTemp.size()<<endl;
					int hhh; cin >> hhh;
				}
			}

			//backup the children before recursive called because the children maybe changed.
			vector<int> vBackUpChildren(vChildren);
			for (vector<int>::iterator citer = vBackUpChildren.begin(); citer != vBackUpChildren.end(); ++citer)
			{
				RoLCA(*citer);
			}
		}
	}

	void setCountDescendant()
	{
		//start from the last deleted vertex, which is the root
		vector<int> ord(iVNum, 0);
		for (int ii = 0; ii < iVNum; ++ii)
		{
			ord[pi[ii]] = ii; //ord[jj] = ii : the vertex ii are the jj=pi[ii]-th deleted vertex
		}
		for (int jj = 0; jj < iVNum; ++jj)
		{
			int vid = ord[jj];


			TreeNode & tn = TG[vid];
			if (tn.children.size() == 0)
				continue;

			int sum = 1;//the node itself
			for (vector<int>::iterator iter = tn.children.begin(); iter != tn.children.end(); ++iter)
			{
				sum += TG[*iter].iCountDescendant;
			}
			tn.iCountDescendant = sum;
		}

		//cout << "root has descendants: " << TG[root].iCountDescendant << " |V|=" << iVNum << endl;
	}

	void resetHeight(int nid)
	{
		//nid = TG[nid].iReal;

		int p = TG[nid].parent;
		if (p != -1)
			TG[nid].height = TG[p].height + 1;
		else
			TG[nid].height = 0;

		for (vector<int>::iterator iter = TG[nid].children.begin(); iter != TG[nid].children.end(); ++iter)
			resetHeight(*iter);
	}

	void writeGoodBranch(FILE * fout)
	{
		fwrite(&Node::s_outputThredhold, sizeof(double), 1, fout);
		int topk = vGoodBranch.size();
		fwrite(&topk, sizeof(int), 1, fout);
		fwrite(vGoodBranch.data(), sizeof(int), topk, fout);
		fwrite(vProbBranch.data(), sizeof(double), topk, fout);
		fwrite(vBetterBranch.data(), sizeof(double), topk, fout);
	}

	void write(string fileOut)
	{
		FILE * fout;
		//fopen_s(&fout, fileOut.c_str(), "wb");
		fout = fopen(fileOut.c_str(), "wb");
		fwrite(&iVNum, sizeof(int), 1, fout);

		writeGoodBranch(fout);
		map<int, double> mBidBetter;
		for(int jj = 0; jj < vGoodBranch.size(); ++jj)
		{
			mBidBetter[vGoodBranch[jj]]=vBetterBranch[jj];
		}
		
		for (int ii = 0; ii < iVNum; ++ii)
			TG[ii].printMe(TG, ii, fout, mBidBetter);


		//write the part for lca
		//fwrite(lca.toRMQ.data(), sizeof(int), iVNum, fout);
		fwrite(lca.toRMQ, sizeof(int), iVNum, fout);

		int iRow = lca.RMQIndex.size();
		fwrite(&iRow, sizeof(int), 1, fout);
		for (int ii = 0; ii < iRow; ++ii)
		{
			int iCol = lca.RMQIndex[ii].size();
			fwrite(&iCol, sizeof(int), 1, fout);
			fwrite(lca.RMQIndex[ii].data(), sizeof(int), iCol, fout);
		}
		fclose(fout);
	}

	void runOneCut2()
	{
		vector<int> vPriority(iVNum, -1);

		int vidLeft[] = {25151, 56600, 56590, 217679, 234240, 233429};
		int vidRight[] = {34160, 56599, 189752, 147077, 225403, 239207};

		int iPriotiry = iVNum;//set to be a large value, delete them at last
		for(int ii=0; ii<6; ++ii)
		{
			vPriority[vidLeft[ii]] = iPriotiry;
		}


		int vidLeft2[] = {59702,59702,62410,62410,64353,64353,64403,64403,64402,64402,64415,64415,64397,64397,64410,64410,64839,64839,64900,64900,65039,65039,65045,65045,65040,65040,65045,65045,65109,65109,65113,65113,65117,65117,65117,65117,65319,65319,65323,65323,65323,65323,65347,65347,72299,72299,72300,72300,74269,74269,62412,62412,64840,64840,77450,77450,63747,63747,78547,78547,78837,78837,78838,78838,62521,62521,85118,85118};
		int vidRight2[] = {59864,59864,62266,62266,64352,64352,64402,64402,64407,64407,64397,64397,64403,64403,64406,64406,64841,64841,64901,64901,62516,62516,65046,65046,65044,65044,65044,65044,65048,65048,65112,65112,65105,65105,65112,65112,65318,65318,65104,65104,65318,65318,65348,65348,72297,72297,72298,72298,72303,72303,76869,76869,64839,64839,62520,62520,65108,65108,78182,78182,62518,62518,77449,77449,62408,62408,90640,90640};
		
		for(int ii=0; ii<68; ++ii)
		{
			vPriority[vidLeft2[ii]] = iPriotiry - 1;
		}

		int vidLeft3[] = {139535,139535,139801,139801,139917,139917,139918,139918,139923,139923,139924,139924,139916,139916,139923,139923,140628,140628,189308,189308,190537,190537,194527,194527,194529,194529,194531,194531,194533,194533,194547,194547,194554,194554,194566,194566,197289,197289,197292,197292,197300,197300,197324,197324,197326,197326,197330,197330,197425,197425,197276,197276,199616,199616,199633,199633,199636,199636,199647,199647,199648,199648};
		int vidRight3[] = {139534,139534,139800,139800,139915,139915,139915,139915,139922,139922,139916,139916,139923,139923,151604,151604,187088,187088,189628,189628,190536,190536,194526,194526,194528,194528,194530,194530,194537,194537,194545,194545,194553,194553,194546,194546,197288,197288,197287,197287,197299,197299,189400,189400,197325,197325,197329,197329,197424,197424,189401,189401,199271,199271,199635,199635,199635,199635,199646,199646,199646,199646};
		
		for(int ii=0; ii<62; ++ii)
		{
			vPriority[vidLeft3[ii]] = iPriotiry - 1;
		}

		treeDecOrdDeg(vPriority);
	}

	bool detect(int eid, Triangle & triangle)
	{
		int start = eh.edges[eid].start;
		int end = eh.edges[eid].end;

		if(vh[start].deleted || vh[end].deleted)
			return false;

		set<int> & edge_list = vh[start].edge_list;
		bool bFound = false;
		for (set<int>::iterator iter = edge_list.begin(); iter != edge_list.end(); ++iter)
		{
			//get id of neighbor
			int vneighbor = eh.edges[*iter].adjVertex(start);

			if (vneighbor > start || vneighbor > end || vh[vneighbor].deleted == true)
				continue;

			//check whether there exist edge (vneighbor, end)
			int thirdEdge = 0;
			if (eh.existEdge(vneighbor, end, thirdEdge))
			{
				if(bFound == false)
				{
					bFound = true;
					triangle.A = start; triangle.B = end; triangle.C = vneighbor;
					triangle.AB = eid; triangle.AC = *iter; triangle.BC = thirdEdge;
				}
				else
				{
					if(vh[vneighbor].edge_list.size() < vh[triangle.C].edge_list.size())
					{
						triangle.C = vneighbor;
						triangle.AC = *iter; triangle.BC = thirdEdge;
					}
				}
			}
		}
		return bFound;
	}

	void handleTriange(Triangle & triangle)
	{
		//remove edges AB,AC,BC
		//add a new vertex D, add edges AD,BD,CD
		LENGTH x = eh.edges[triangle.AB].length, y = eh.edges[triangle.AC].length, z = eh.edges[triangle.BC].length;
		LENGTH a = (x+y-z)/2, b = (x+z-y)/2, c = (y+z-x)/2;

		vh[triangle.A].deleteANeighborEdge(triangle.AB);vh[triangle.A].deleteANeighborEdge(triangle.AC);
		vh[triangle.B].deleteANeighborEdge(triangle.AB);vh[triangle.B].deleteANeighborEdge(triangle.BC);
		vh[triangle.C].deleteANeighborEdge(triangle.AC);vh[triangle.C].deleteANeighborEdge(triangle.BC);

		eh.removeAnEdge(triangle.AB); eh.removeAnEdge(triangle.AC); eh.removeAnEdge(triangle.BC);

		Vertex Dtemp;
		Dtemp.id=vh.size();
		
		VertexInH D(Dtemp);

		Edge AD, BD, CD;
		AD.start = triangle.A; AD.end = D.id; AD.length = a;
		BD.start = triangle.B; BD.end = D.id; BD.length = b;
		CD.start = triangle.C; CD.end = D.id; CD.length = c;

		int eidAD = eh.edges.size(), eidBD = eidAD+1, eidCD = eidAD+2;
		eh.insertAnEdge(AD); eh.insertAnEdge(BD); eh.insertAnEdge(CD);

		D.edge_list.insert(eidAD); vh[triangle.A].insertANeighborEdge(eidAD);
		D.edge_list.insert(eidBD); vh[triangle.B].insertANeighborEdge(eidBD);
		D.edge_list.insert(eidCD); vh[triangle.C].insertANeighborEdge(eidCD);


		vh.push_back(D);
	}

	void dealTriangle()
	{
		const int iEdgeNum = eh.edges.size();
		int iCountTriangle = 0;
		for(int eid = 0; eid < iEdgeNum; ++eid)
		{
			if(eh.edges[eid].mark == -1)
				continue;
			Triangle triangle;
			if(detect(eid, triangle))
			{
				++iCountTriangle;
				handleTriange(triangle);
			}
		}
		cout<<"Num of Triangles = "<<iCountTriangle<<endl;
		iVNum = vh.size();
	}

	void showTreeStructure(int iRoot, int iCurrent, int & iCount)
	{
		int iNumOfChildren = TG[iCurrent].children.size();

		if(iNumOfChildren == 0)
		{
			++iCount;
			if(iCount > 10)
				cout<<"\t"<<iRoot<<","<<iCount<<endl;
		}
		else if(iNumOfChildren == 1)
		{
			++iCount;
			showTreeStructure(iRoot, TG[iCurrent].children[0], iCount);
		}
		else
		{
			++iCount;
			if(iCount > 10)
				cout<<"\t"<<iRoot<<","<<iCount<<endl;
			for(int ii = 0; ii < iNumOfChildren; ++ii)
			{
				int count = 0;
				showTreeStructure(TG[iCurrent].children[ii], TG[iCurrent].children[ii], count);
			}
		}
	}

	double avgHeight()
	{
		double sum = 0;
		for(int ii = 0; ii < iVNum; ++ii)
			sum += TG[ii].height;
		cout<<"avgHeight="<<sum/iVNum<<'\t';
		return sum/iVNum;
	}
	double avgWidth()
	{
		setCountDescendant();
		double width = 0.0;
		double count = 0.0;
		int maxW = 0;
		for(int ii = 0; ii < iVNum; ++ii)
		{
			const int iSize = TG[ii].children.size();
			if(iSize < 2)
				continue;
			double iPairs = 0;
			const vector<int> & children = TG[ii].children;
			int sum = TG[ii].iCountDescendant - 1;
			for(int jj = 0; jj < iSize; ++jj)
			{
				int ds = TG[children[jj]].iCountDescendant;
				sum -= ds;
				iPairs += 1.0*ds*sum;
			}
			assert(sum == 0);

			width += iPairs * TG[ii].pos.size();

			count += iPairs;

			if(TG[ii].pos.size()>maxW)
				maxW = TG[ii].pos.size();
		}
		width = width/count;
		cout<< "avgWidth="<<width << '\t'<<"maxW="<<maxW<<endl;
		mResult['a'] = width; mResult['w'] = maxW;
		return width;
	}
	void setBranchHeight(int bid, int bheight)
	{
		const int iSize = TG[bid].children.size();
		if(iSize > 1)
		{
			TG[bid].branch_height = bheight;
			for(vector<int>::iterator iter = TG[bid].children.begin(); iter != TG[bid].children.end(); ++iter)
			{
				setBranchHeight(*iter, bheight+1);
			}
		}
		else if(iSize == 1)
		{
			setBranchHeight(TG[bid].children[0], bheight);
		}
	}

	void showWidthTendency()
	{
		int iThreshold = 0.1*avgHeight();
		int iCount = 0;
		int iTop = root;
		while(TG[iTop].children.size()==1)
			iTop = TG[iTop].children[0];
		for(int ii = 0; ii < iVNum; ++ii)
		{
			if(TG[ii].children.size() != 0)
				continue;
			int count = 0;
			int current = ii;
			int w = TG[current].neighbors.size();
			int parent = TG[current].parent;

			while(parent != iTop)
			{
				int w1 = TG[parent].neighbors.size();
				if(w1 > w)
					++count;
				w = w1;
				current = parent;
				parent = TG[current].parent;
			}
			if(count > iThreshold)
				++iCount;
		}
		cout<<"There are "<<iCount<<" branches where more than "<<iThreshold<<" nodes decrease the width"<<endl;
	}

	void handleDescendants(const int bid, int current, const vector<int>& vCuts, const vector<int>& vCutsHeight)
	{
		//for each desecendant v, prune the pos by doing checking on vCuts.
		//d(v, ci) >= d(v, cj) + d(cj, ci)
		const int iSize = vCutsHeight.size();
		
		for(vector<int>::iterator iter = TG[current].children.begin(); iter != TG[current].children.end(); ++iter)
		{
			const int v = *iter;
			const LENGTH * vToAncestor = TG[v].vToAncestor.data();
			for(int ii = 0; ii < iSize; ++ii)
			{
				bool unpruned = true;
				const LENGTH dvci = vToAncestor[vCutsHeight[ii]];
				for(int jj = 0; jj < iSize; ++jj)
				{
					if(ii == jj)
						continue;
					LENGTH dvcj = vToAncestor[vCutsHeight[jj]];
					if(dvci <= dvcj)
						continue;
					LENGTH dcjci = 0;
					if(vCutsHeight[ii]<vCutsHeight[jj])
						dcjci = TG[vCuts[jj]].vToAncestor[vCutsHeight[ii]];
					else
						dcjci = TG[vCuts[ii]].vToAncestor[vCutsHeight[jj]];
					if(dvci >= dvcj + dcjci)
					{
						unpruned = false;
						break;
					}
				}
				if(unpruned)
					TG[v].mPos[bid].push_back(vCutsHeight[ii]);

			}
			handleDescendants(bid, *iter, vCuts, vCutsHeight);
		}
	}

	void PLLABranchNode(int bid)
	{
		const vector<int> & vNeighbors=TG[bid].neighbors;
		const int iNum = vNeighbors.size() ;
		const vector<int> & vPos=TG[bid].pos;
		const int iNum2 = vPos.size();//pos is the height of ancestors, also the cut
		
		vector<bool> vUsed(maxHeight, false);
		 for(int jj = 0; jj < iNum2; ++jj)
		 	vUsed[vPos[jj]]=true;

		vector<int> vCuts;
		vector<int> vCutsHeight;
		for(int ii = 0; ii < iNum; ++ii)
		{
			int height = TG[vNeighbors[ii]].height;
			if(vUsed[height])
			{
				vCuts.push_back(vNeighbors[ii]);
				vCutsHeight.push_back(height);
			}
		}

		handleDescendants(bid, bid, vCuts, vCutsHeight);
	}
	void buildAllBranches()
	{
		for (int jj = 0; jj < vGoodBranch.size(); ++jj)
		{
			int vid = vGoodBranch[jj];
			PLLABranchNode(vid);
		}
	}

	void getProbForBranchNode(int topk = 20)
	{
		const int iTotal = g_Samples;
		vector<int> vCountLca(iVNum, 0);
		int sum = 0;
		for(int ii = 0; ii < iTotal; ++ii)
		{
			int v1 = rand() / (double)RAND_MAX*iVNum;
			int v2 = rand() / (double)RAND_MAX*iVNum;
			int vc = lca.LCAQuery(v1,v2);
			if(vc == v1 || vc == v2)
				continue;
			vCountLca[vc] +=1;
			sum += 1;
		}
		cout<<sum<<" pairs need to find lca among "<<iTotal<<" pairs"<<endl;

		MinHeap mh;
		mh.initial(topk, iVNum);

		for(int ii = 0; ii < topk; ++ii)
			mh.Push(vCountLca[ii], ii);

		int iTop = mh.top();
		
		for(int ii = topk; ii < iVNum; ++ii)
		{
			if(vCountLca[ii] > iTop)
			{
				int temp;
				iTop = mh.Pop(temp);
				mh.Push(vCountLca[ii], ii);
			}
		}

		double psum = 0;
		vGoodBranch.resize(topk);
		vProbBranch.resize(topk);
		for(int ii = 0; ii < topk; ++ii)
		{
			int vid, count;
			count = mh.Pop(vid);
			vGoodBranch[topk-1-ii] = vid;
			vProbBranch[topk-1-ii] = count*1.0/sum;
			assert(vCountLca[vid]==count);
			cout<<"vid="<<vid<<"\t prob="<<count*1.0/sum<<"\t branch_height="<<TG[vid].branch_height<<endl;
			psum += count*1.0/sum;
		}
		cout << "The top "<<topk<<" branch nodes covers "<<psum*100<<"%."<<endl;
		//for(int ii = 0; ii < topk; ++ii)
			//cout<<"vid="<<vGoodBranch[ii]<<endl;
	}

	double resizeBetterBranch(int rsize = 5)
	{
		//pick top-k from betterbranch
		MinHeap mh;
		const int topk = vGoodBranch.size();
		mh.initial(topk, topk);

		for(int ii = 0; ii < rsize; ++ii)
			mh.Push((int)(vBetterBranch[ii]*10000), ii);

		int iTop = mh.top();
		
		for(int ii = rsize; ii < topk; ++ii)
		{
			if((int)(vBetterBranch[ii]*10000) > iTop)
			{
				int temp;
				iTop = mh.Pop(temp);
				mh.Push((int)(vBetterBranch[ii]*10000), ii);
			}
		}

		vector<int> vBackupGB(vGoodBranch);
		vector<double> vBackupBB(vBetterBranch);
		vector<double> vBackupPB(vProbBranch);


		vGoodBranch.resize(rsize);
		vProbBranch.resize(rsize);
		double decrease = 0;
		for(int ii = 0; ii < rsize; ++ii)
		{
			int index, count;
			count = mh.Pop(index);
			//cout<<count<<'\t';
			int jj = rsize-1-ii;
			vGoodBranch[jj] = vBackupGB[index];
			vBetterBranch[jj] = vBackupBB[index];
			vProbBranch[jj] = vBackupPB[index];
			decrease += vBetterBranch[jj];
			//cout << vBetterBranch[jj]<<'\t';
		}
		return decrease;
	}


	void tellBetterBranch()
	{
		//compare the decrease on average width
		const int iTotal = g_Samples;
		const int iSize = vGoodBranch.size();
		vBetterBranch.resize(iSize, 0);
		vector<double> vSumWidth(iSize, 0);
		vector<int> vCount(iSize, 0);
		int sum = 0;
		double avgpw = 0;
		double avgpw2=0;
		int iPairs = 0;
		for(int ii = 0; ii < iTotal; ++ii)
		{
			int v1 = rand() / (double)RAND_MAX*iVNum;
			int v2 = rand() / (double)RAND_MAX*iVNum;
			int vc = lca.LCAQuery(v1,v2);
			if(vc == v1 || vc == v2)
				continue;

			++iPairs;
			avgpw2 += TG[vc].pos.size();
				
			bool bPrune = false;
			for(int jj = 0; jj < iSize; ++jj)
			{
				if(vc == vGoodBranch[jj])
				{
					bPrune = true;
					++vCount[jj];
					if((TG[v1].mPos)[vc].size() < (TG[v2].mPos)[vc].size())
					{
						vSumWidth[jj] += (TG[v1].mPos)[vc].size();
						avgpw += (TG[v1].mPos)[vc].size();
					}
					else
					{
						vSumWidth[jj] += (TG[v2].mPos)[vc].size();
						avgpw += (TG[v2].mPos)[vc].size();
					}
					break;
				}
			}
			if(!bPrune)
			{
				avgpw += TG[vc].pos.size();
			}
		}
		for(int jj = 0; jj < iSize; ++jj)
		{
			if(vCount[jj] == 0)
				continue;
			int temp = vSumWidth[jj] / vCount[jj];

			int width = TG[vGoodBranch[jj]].pos.size();
			vBetterBranch[jj] = (width - temp)*vProbBranch[jj] ;
		}
		for(int ii = 0; ii < iSize; ++ii)
			cout<<"vid="<<vGoodBranch[ii]<<"\tprob="<<vProbBranch[ii]<<"\t width_discount="<<vBetterBranch[ii]<<"\t branch_height="<<endl;

		cout <<"avgpw="<<avgpw/iPairs<<'\t'<<"avgpw2="<<avgpw2/iPairs<<'\t';

		double decrease = resizeBetterBranch();

		cout <<"avgpw-topk = " << avgpw2/iPairs - decrease<<endl;
		mResult['a'] = avgpw2/iPairs - decrease;
	}

	void setPriority(vector<int> & vPriority)
	{
		const int iTotal = g_Samples;
		vPriority.resize(iVNum, -1);


		for(int ii = 0; ii < iTotal; ++ii)
		{
			int v1 = rand() / (double)RAND_MAX*iVNum;
			int v2 = rand() / (double)RAND_MAX*iVNum;
			int vc = lca.LCAQuery(v1,v2);
			if(vc == v1 || vc == v2)
				continue;

			int ps = TG[vc].pos.size();
			int * p2 = TG[vc].pos.data();
			LENGTH *dx = TG[v1].vToAncestor.data(), *dy = TG[v2].vToAncestor.data();

			LENGTH res = (int)INFINITY;
			//int index = 0;
			vector<int> vIndex;
			for (int i = 0; i < ps; i++){
				LENGTH tmp = dx[p2[i]] + dy[p2[i]];
				if (res > tmp)
				{
					vIndex.resize(0);
					vIndex.push_back(i);
					res = tmp;
				}
				else if(res == tmp)
				{
					vIndex.push_back(i);
				}
			}

			if(vIndex.size() > 0)
			{
				for(int jj = 0; jj < vIndex.size(); ++jj)
				{
					int index = vIndex[jj];
					if(index < ps - 1)
					{
						int vid = TG[vc].neighbors[index];
						++vPriority[vid];
					}

				}
				
			}
			else
				++vPriority[vc];
		}
	}

	inline double getLearningRate(double oldw, double neww)
	{
		//if oldw < neww, set the rate small
		//else set the rate large
		double lr;
		if(oldw < neww)
			lr = 0.05;
		else
			lr = 0.95;
		return lr;
	}

	void updatePriority(vector<int> & vopt, double & optw, vector<int> & vcurr, double currw, const vector<int> & vnew, double & learning_rate)
	{
		cout<<" learning_rate="<<learning_rate<<"\t optw="<<optw<<"\t currw="<<currw<<'\t';
		if(optw >= currw)
		{
			optw = currw;
			//learning_rate = (learning_rate >= 0.2 ? learning_rate*2 : 0.2) ;
			for(int ii = 0; ii < iVNum; ++ii)
			{
				vopt[ii] = vcurr[ii];
				vcurr[ii] = vnew[ii];
				//vcurr[ii] = vcurr[ii] + learning_rate * (vnew[ii] - vcurr[ii]);
			}
		}
		else
		{
			learning_rate = 1;
			for(int ii = 0; ii < iVNum; ++ii)
			{
				vcurr[ii] = vnew[ii];
				//vcurr[ii] = vopt[ii] + learning_rate * (vnew[ii] - vopt[ii]);
			}
		}
	}
	void boost(int iBoostTimes)
	{
		initial();
		treeDecOnlyDegree();
		maketree();

		calDisAndPos(root);
	
		if(g_bRoLCA)
		{
			//cout<<endl<<"before RoLCA:\t";
					//avgHeight();
		//avgWidth();
			RoLCA(root);
			//cout<<endl<<"after RoLCA:\t";
		}

		avgHeight();
		double opt_width = avgWidth();
		double curr_width = opt_width;
		vector<int> vOptPriority(iVNum, -1);
		vector<int> vCurrPriority(iVNum, -1);
		double learning_rate = 0.1;

		setBranchHeight(root, 0);

		lca.makeRMQ(root);

		//getProbForBranchNode(bheight_control);
		//buildAllBranches();
		//tellBetterBranch();


		for(int ii = 0; ii < iBoostTimes; ++ii)
		{
			vector<int> vPriority(iVNum, -1);
			setPriority(vPriority);/*
			{
				RoLCA(root);
				cout<<endl<<"after RoLCA:\t";
				avgWidth();
			}*/
			//learning_rate = 1.0/(1+ii);
			updatePriority(vOptPriority, opt_width, vCurrPriority, curr_width, vPriority, learning_rate);
			initial();
			treeDecDegOrd(vCurrPriority);

			maketree();

			calDisAndPos(root);
		
			if(g_bRoLCA)
			{
				//cout<<endl<<"before RoLCA:\t";
									//avgHeight();
		//avgWidth();
				RoLCA(root);
				//cout<<endl<<"after RoLCA:\t";
			}

			avgHeight();
			curr_width = avgWidth();
			//learning_rate = getLearningRate(old_width, new_width);
			setBranchHeight(root, 0);
			lca.makeRMQ(root);

		}
		if(curr_width > opt_width)
		{
			initial();
			treeDecDegOrd(vOptPriority);
			maketree();
		}
	}

	void showTraditionalWidth()
	{
		double avgw = 0;
		double avgb = 0;
		int countb = 0;
		double maxbHeight = 0;
		for(int ii = 0; ii < iVNum; ++ii)
		{
			avgw += TG[ii].pos.size();
			if(TG[ii].children.size()>1)
			{
				avgb += TG[ii].pos.size();
				++countb;
				if(TG[ii].branch_height > maxbHeight)
					maxbHeight = TG[ii].branch_height;
			}
		}
		cout<<"\nTraditional: avgw = "<<avgw/iVNum <<" avgb = "<<avgb/countb <<" maxbHeight = "<<maxbHeight<<endl;
	}

	void output(string dataset, ofstream & fout)
	{
		fout << dataset << '\t';
		fout << mResult['o'] <<'\t' << mResult['r'] << '\t' << mResult['p'] << '\t';
		fout << mResult['a'] << '\t' << mResult['w'] << '\t' << mResult['t'] << '\t' << mResult['m'] << '\t';
		fout << mResult['L'] << '\t' << mResult['R'] << '\t' << mResult['P'] << endl;
	}

	vector<int> vMarkInHeap;
	MinHeap2 minheap;
	double singleSource(int vid, int target, int mark)
	{		
		KEY_TYPE key(0,0,0);
		minheap.Push(key, vid);
		vMarkInHeap[vid] = mark;

		while(minheap.count)
		{
			int vtop;
			KEY_TYPE currentKey = minheap.Pop(vtop);
			if(vtop == target)
				return currentKey.key1;

			//for the neighbor of the vertex
			for (set<int>::iterator jj = vh[vtop].edge_list.begin(); jj != vh[vtop].edge_list.end(); ++jj)
			{
				int eid = *jj;
				int vneighbor = eh.edges[eid].adjVertex(vtop);
				int pos = minheap.index[vneighbor];
				if(pos == minheap.NO)
					continue;
				if (mark == vMarkInHeap[vneighbor])
				{				
					//2. out of heap
					assert(pos < minheap.count);					
					//1. still in minheap
					KEY_TYPE newkey(currentKey.key1+eh.edges[eid].length, 0, 0);
					minheap.UpdateKey(pos, newkey);
				}
				else
				{
					KEY_TYPE newkey(currentKey.key1+eh.edges[eid].length, 0, 0);
					minheap.Push(newkey, vneighbor);
					vMarkInHeap[vneighbor]=mark;				
				}
			}
		}
	}

	void testD(int iCount)
	{
		initial();
		vector<int> v1(iCount);
		vector<int> v2(iCount);
		srand(time(NULL));
		for(int ii = 0; ii < iCount; ++ii)
		{
			v1[ii] = rand()*1.0/RAND_MAX * iVNum;
			v1[ii] = rand()*1.0/RAND_MAX * iVNum;
		}
		double dSum = 0;
		int timestamp = 0;
		vMarkInHeap.resize(iVNum,-1);
		minheap.initial(iVNum);
		clock_t start = clock();
		for(int ii = 0; ii < iCount; ++ii)
		{
			minheap.Reset();
			++timestamp;
			dSum += singleSource(v1[ii], v2[ii], timestamp);
		}
		cout<<endl;
		cout<< (clock()-start)*1.0/iCount << '\t' << dSum/iCount<<endl;
	}

};
double Node::s_outputThredhold = 0.5;

int main(int argc, char* argv[])
{
	string dataPath("../data/");
	string dataset("NY");

	int method = 3;//

	bool bRoLCA = true;

	int bheight_control = 50;

	int ii = 4;

	Node::s_outputThredhold = 0;

	if(argc>ii)
		bheight_control = atoi(argv[ii]);

	--ii;
	if(argc > ii)
		bRoLCA = atoi(argv[ii]);

	--ii;
	if (argc>ii)
		method = atoi(argv[ii]);

	--ii;
	if (argc>ii)
		dataset = argv[ii];

	cout << dataset << '\t' << method << '\t';

	string dataname("USA-road-d." + dataset + ".gr");
	string coordinate("USA-road-d." + dataset + ".co");

	string filename(dataPath + dataname);
	timer t0;
	t0.restart();
	clock_t start = clock();

	ReadFile rf(filename);

	rf.readRoadNet(filename);

	char indexName[100];
	sprintf(indexName, "%s_%d-%d-%d.pindex", dataset.c_str(), method, bRoLCA, bheight_control);

	if(method == COOR)
		rf.readCoordinate(dataPath+coordinate);

	//rf.edgeList.component(47437, rf.vertexList);

	Tree_Decomposition td(rf);
	//td.testD(method);
	//return 0;
	int btimes = 3;
	if(dataset == "CAL")
		btimes = 1;
	if(dataset == "FLA")
		btimes = 4;
	if(dataset == "E")
		btimes = 2;

	td.g_bRoLCA = bRoLCA;
	clock_t startLearn = clock();
	if(method == LB)
	{
		td.g_Samples = 10000;
		//cout<<td.g_Samples<<endl;
		td.boost(5);
		//return 0;
	}
	else
		td.run(method);
	td.mResult['L']=(double)(clock() - startLearn) / CLOCKS_PER_SEC;//time for learning

	timer t1;
	t1.restart();
	td.calDisAndPos(td.root);
	cout<<"caldist takes "<<t1.getTime()<<endl;
	
	clock_t startRoLCA = clock();
	timer t2;
	t2.restart();
	if(bRoLCA)
		td.RoLCA(td.root);
	td.mResult['R']= (double)(clock() - startRoLCA) / CLOCKS_PER_SEC;//time for RoLCA
	cout<<"rolca takes"<<t2.getTime()<<endl;
	
	td.avgHeight();
	td.avgWidth();

	td.setBranchHeight(td.root, 0);

	td.lca.makeRMQ(td.root);
	clock_t startProject = clock_t();
	timer t3;
	t3.restart();
	if(bheight_control > 0 && bheight_control < 100){
		td.g_Samples = 1000000;
		td.getProbForBranchNode(bheight_control);
		td.buildAllBranches();
		td.tellBetterBranch();
	}
	cout<<"project takes "<<t3.getTime()<<endl;
	td.mResult['P']= (double)(clock() - startProject) / CLOCKS_PER_SEC;//time for project
	
	//td.showTraditionalWidth();
	clock_t startwrite=clock();
	td.write(indexName);
	cout <<"write needs="<< (double)(clock() - startwrite) / CLOCKS_PER_SEC << '\t';
	
	td.mResult['t'] = (double)(clock() - start) / CLOCKS_PER_SEC;
	cout << (double)(clock() - start) / CLOCKS_PER_SEC << '\t';
	cout << "overall takes "<<t0.getTime()<<endl;
	string indexFile(indexName);

	td.mResult['m'] = file_size(indexFile.c_str());
	//cout << indexSize/1000.0/1000.0 << '\t';
	cout << endl;
	//return 0;
	td.mResult['o'] = method;
	td.mResult['r'] = bRoLCA;
	td.mResult['p'] = bheight_control;
	ofstream fout;
	fout.open("records.txt", ios::app);
	td.output(dataset,fout);
	fout.close();

	return 0;
}
