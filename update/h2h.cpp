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
	string dataset;
	int g_bRoLCA;
	Tree_Decomposition(ReadFile & rf,string _dataset="") :pvl(&rf.vertexList), pel(&rf.edgeList),g_bRoLCA(false),dataset(_dataset)
	{
	}
	~Tree_Decomposition()
	{
		//cout<<pvl->size()<<endl;
		//cout<<pel->edges.size()<<endl;
	}


	void run(int method)
	{
		initial();


		switch (method)
		{
		case DegreeOnly:
			treeDecOnlyDegree();
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

	vector<VertexInH> vh;
	int iVNum, iENum, iENum0;
	EdgeList eh;
	vector<TreeNode> TG;
	vector<int> pi;
	int order;
	int root;

	LCA lca;

	vector<set<TRIPLE> > TabE;
	vector<LENGTH> TabPhi;

	void initial()
	{
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

		iENum0=iENum=eh.edges.size();
		TabE.resize(iENum);
		TabPhi.resize(iENum);
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

			vh[vid].removeMe(vh, eh, tn, false, TabE, TabPhi);

			//update the degree of its neighbors
			for (int jj = 0; jj < tn.neighbors.size(); ++jj)
			{
				int vneighbor = tn.neighbors[jj];
				int degree = vh[vneighbor].edge_list.size();

				int pos = minheap.index[vneighbor];
				minheap.UpdateKey(pos, degree);
			}

		}
		cout<<endl<<TabPhi.size()<<'\t'<<TabE.size()*1.0/iENum0<<endl;
		double sum = 0;
		for(int ii=0; ii < TabE.size(); ++ii)
		{
			int isize = TabE[ii].size();
			sum += isize;
		}
		cout<<"avgTabE="<<sum/TabE.size()<<endl;
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

	void maketree()
	{
		iENum = eh.edges.size();
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
		int maxHeight = 0;
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

	void write(string fileOut)
	{
		FILE * fout;
		//fopen_s(&fout, fileOut.c_str(), "wb");
		fout = fopen(fileOut.c_str(), "wb");
		fwrite(&iVNum, sizeof(int), 1, fout);

		for (int ii = 0; ii < iVNum; ++ii)
			TG[ii].printMe(TG, ii, fout);


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

	void writeTable(string fileOut)
	{
		FILE * fout;
		fout = fopen(fileOut.c_str(), "wb");
		fwrite(TabPhi.data(), sizeof(LENGTH), iENum, fout);

		for(int ii = 0; ii < iENum; ++ii)
		{
			set<TRIPLE> & triple=TabE[ii];
			int iSize=triple.size();
			fwrite(&iSize, sizeof(int), 1, fout);
			for(set<TRIPLE>::iterator iter=triple.begin(); iter!=triple.end(); ++iter)
			{
				fwrite(&(iter->vid), sizeof(int), 1, fout);
				fwrite(&(iter->eid1), sizeof(int), 1, fout);
				fwrite(&(iter->eid2), sizeof(int), 1, fout);
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

	void output(string dataset, ofstream & fout)
	{
		fout << dataset << '\t';
		fout << mResult['o'] <<'\t' << mResult['r'] << '\t' << mResult['p'] << '\t';
		fout << mResult['a'] << '\t' << mResult['w'] << '\t' << mResult['t'] << '\t' << mResult['m'] << endl;
	}

	void updateDis(int iRoot, int timestamp)
	{
		//should not be called before the height is set correctly
		//calculate the distance to each ancestors

		TreeNode & tn = TG[iRoot];
		if(timestamp == tn.iTimeStamp)
			return;
		tn.iTimeStamp = timestamp;

		tn.vToAncestor.resize(tn.height + 1);
		tn.vToAncestor[tn.height] = 0;
		vector<int> & vNeighbor = tn.neighbors;
		int ancestor = tn.parent;
		const int iN = vNeighbor.size();
		vector<LENGTH> vphi(iN);
		for(int ii = 0; ii < iN; ++ii)
		{
			int nvid = vNeighbor[ii];
			int eid;
			eh.existEdge(iRoot, nvid, eid);
			vphi[ii] = TabPhi[eid];
		}
		while (ancestor != -1)
		{
			//calculate the distance to the ancestor
			//enumerate all pathes through the neighbors
			LENGTH minDis = (int)INFINITY;
			int posOfAncestor = TG[ancestor].height;

			for (int ii = 0; ii < iN; ++ii)
			{
				int nvid = vNeighbor[ii];
				LENGTH temp = vphi[ii];
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

		//recursive on children
		for (int ii = 0; ii < tn.children.size(); ++ii)
			updateDis(tn.children[ii], timestamp);
	}

	MinHeap2 QQ;
	double dCaldistTime;
	void update(int eid, int u, int v, LENGTH ww)
	{
		if(pi[u]>pi[v])
		{
			update(eid, v, u, ww);
			return;
		}
		LENGTH newphi=(LENGTH)INFINITY;
		if(ww < TabPhi[eid])
			newphi = ww;
		else
		{
			const set<TRIPLE> & triple = TabE[eid];
			
			for(set<TRIPLE>::iterator iter=triple.begin(); iter!=triple.end(); ++iter)
			{
				LENGTH temp = TabPhi[iter->eid1]+TabPhi[iter->eid2];
				if(temp < newphi)
					newphi = temp;
			}
			newphi=min(ww, newphi);
			if(newphi == TabPhi[eid])
			{
				eh.edges[eid].length = ww;
				return;
			}
		}
		int r=u;
		KEY_TYPE::s_pi = &pi;		
		
		QQ.Reset();
		KEY_TYPE key(u,v, newphi);
		if(!eh.existEdge(u,v,eid))
			assert(false);
		QQ.Push(key, eid);
		TabPhi[eid] = newphi;
		eh.edges[eid].length = ww;
		while(QQ.count)
		{
			vector<pair<KEY_TYPE, int> > vQ;
			vQ.push_back(QQ.minheap[0]);
			map<int, LENGTH> mOldPhi;
			int eid0;
			KEY_TYPE currentKey = QQ.Pop(eid0);
			mOldPhi[eid0] = TabPhi[eid0];
			TabPhi[eid0]=currentKey.key3;
			int v1 = currentKey.key1;
			r = v1;
			
			while(QQ.count)
			{
				if(QQ.minheap[0].first.key1 == v1)
				{
					vQ.push_back(QQ.minheap[0]);
					int temp;
					KEY_TYPE ktemp=QQ.Pop(temp);
					mOldPhi[temp] = TabPhi[temp];
					TabPhi[temp] = ktemp.key3;
				}
				else
					break;
			}
			int iSizevQ=vQ.size();
			TreeNode & tn=TG[v1];
			vector<int> & vNeighbor = tn.neighbors;
			map<int, bool> mHandledEdge;
			map<int, int> v2e;
			for (int ii = 0; ii < vNeighbor.size(); ++ii)
			{
				int c = vNeighbor[ii];
				int e1=-1;
				if(!eh.existEdge(v1, c, e1))
					assert(false);
				v2e[c]=e1;
			}

			for (int ii = 0; ii < vNeighbor.size(); ++ii)
			{
				int c = vNeighbor[ii];
				int e1=v2e[c];
				LENGTH phi1_old;
				if(mOldPhi.find(e1)==mOldPhi.end())
					phi1_old=TabPhi[e1];
				else
					phi1_old=mOldPhi[e1];

				for(int jj = 0; jj < iSizevQ; ++jj)
				{
					int y=vQ[jj].first.key2;
					if(c==y)
						continue;
					int e2=v2e[y];
					int e=-1;
					if(!eh.existEdge(c, y, e))
						assert(false);
					if(mHandledEdge.find(e)!=mHandledEdge.end())
						continue;
					mHandledEdge[e]=true;
					LENGTH phi2_old = mOldPhi[e2];

					//should consider whether e is in Q or not
					//if not, if phi remain, ignore, else update.
					//if yes, if phi increase, or remian, ignore; else update phi.
					LENGTH bk=TabPhi[e];
					bool eNotInQ=QQ.notIn(e);
					//assert(Q.index[e]!=Q.NO);//e should not be pop out before

					
					LENGTH newphie=TabPhi[e1]+TabPhi[e2];
					//if the current wedge gives the minimal, then read historical records for update
					//if not, compare the wedge and the old phi
					if(phi1_old+phi2_old == TabPhi[e])
					{
						const set<TRIPLE> & triple = TabE[e];

						for(set<TRIPLE>::iterator iter=triple.begin(); iter!=triple.end(); ++iter)
						{
							if(pi[iter->vid] > pi[v1])
								continue;
							LENGTH temp = TabPhi[iter->eid1] + TabPhi[iter->eid2];
							if(temp < newphie)
								newphie = temp;
						}
					}
					else
						newphie=min(TabPhi[e], newphie);
					if(e < iENum0)
						newphie = min(eh.edges[e].length, newphie);
					if(bk == newphie)
						continue;
					if(eNotInQ)
					{
						TabPhi[e] = newphie;
						if(pi[c]<pi[y])
						{
							KEY_TYPE ktemp(c,y,TabPhi[e]);
							QQ.Push(ktemp, e);
							//if(!eh.existEdge(c, y, e))
						//assert(false);
						}
						else
						{
							KEY_TYPE ktemp(y,c,TabPhi[e]);
							QQ.Push(ktemp, e);
							//if(!eh.existEdge(c, y, e))
						//assert(false);
						}
					}
					else if(bk > newphie)
					{
						TabPhi[e] = newphie;
						//assert(Q.minheap[Q.index[e]].second == e);
						QQ.minheap[QQ.index[e]].first.key3 = newphie;
					}
				}
			}
		}
		static int timestamp=0;
		clock_t start = clock();
		updateDis(r, timestamp);
		dCaldistTime += clock()-start;
		++timestamp;
	}
	void update2(int eid, int u, int v, LENGTH ww)
	{
		if(pi[u]>pi[v])
		{
			update(eid, v, u, ww);
			return;
		}
		LENGTH newphi=(LENGTH)INFINITY;
		if(ww < TabPhi[eid])
			newphi = ww;
		else
		{
			const set<TRIPLE> & triple = TabE[eid];
			
			for(set<TRIPLE>::iterator iter=triple.begin(); iter!=triple.end(); ++iter)
			{
				LENGTH temp = TabPhi[iter->eid1]+TabPhi[iter->eid2];
				if(temp < newphi)
					newphi = temp;
			}
			newphi=min(ww, newphi);
			if(newphi == TabPhi[eid])
			{
				eh.edges[eid].length = ww;
				return;
			}
		}
		int r=u;
		KEY_TYPE::s_pi = &pi;		
		
		QQ.Reset();
		KEY_TYPE key(u,v, newphi);
		if(!eh.existEdge(u,v,eid))
			assert(false);
		QQ.Push(key, eid);
		map<int, LENGTH> phibk;
		phibk[eid] = TabPhi[eid];
		TabPhi[eid] = newphi;
		eh.edges[eid].length = ww;
		while(QQ.count)
		{
			vector<pair<KEY_TYPE, int> > vQ;
			vQ.push_back(QQ.minheap[0]);

			int eid0;
			KEY_TYPE currentKey = QQ.Pop(eid0);
			phibk[eid0] = TabPhi[eid0];
			TabPhi[eid0]=currentKey.key3;

			int v1 = currentKey.key1;
			r = v1;
		
			while(QQ.count)
			{
				if(QQ.minheap[0].first.key1 == v1)
				{
					vQ.push_back(QQ.minheap[0]);
					int temp;
					KEY_TYPE ktemp=QQ.Pop(temp);
					phibk[temp] = TabPhi[temp];
					TabPhi[temp] = ktemp.key3;
				}
				else
					break;
			}
			int iSizevQ=vQ.size();
			TreeNode & tn=TG[v1];
			vector<int> & vNeighbor = tn.neighbors;
			map<int, int> vNeighborEdge;
			map<int, LENGTH> vTmpPhi, vBkPhi;
			for(int ii = 0; ii < vNeighbor.size(); ++ii)
			{
				int c = vNeighbor[ii];
				int e1 = -1;
				if(!eh.existEdge(v1,c,e1))
					assert(false);
				vNeighborEdge[c] = e1;
				vTmpPhi[c] = TabPhi[e1];
				if(phibk.find(e1)!=phibk.end())
					vBkPhi[c] = phibk[e1];
				else
					vBkPhi[c] = vTmpPhi[c];
			}
			for (int ii = 0; ii < vNeighbor.size(); ++ii)
			{
				int c = vNeighbor[ii];
				//int e1=-1;
				//if(!eh.existEdge(v1, c, e1))
					//assert(false);
				LENGTH phi1 = vTmpPhi[c];
				LENGTH phi1_0 = vBkPhi[c];

				for(int jj = 0; jj < iSizevQ; ++jj)
				{
					int y=vQ[jj].first.key2;
					if(c==y)
						continue;
					//int e2=-1;
					//if(!eh.existEdge(v1, y, e2))
					//	assert(false);
					LENGTH phi2 = vTmpPhi[y];
					LENGTH phi2_0 = vBkPhi[y];
					int e=-1;
					if(!eh.existEdge(c, y, e))
						assert(false);
					//should consider whether e is in Q or not
					//if not, if phi remain, ignore, else update.
					//if yes, if phi increase, or remian, ignore; else update phi.
					LENGTH bk=TabPhi[e];
					LENGTH newphie=(LENGTH)INFINITY;
					if(bk < phi1_0+phi2_0)
						newphie = min(bk, phi1+phi2);
					else
					{
						if(phi1+phi2<=bk)
							newphie = phi1+phi2;
						else
						{
							const set<TRIPLE> & triple = TabE[e];
							for(set<TRIPLE>::iterator iter=triple.begin(); iter!=triple.end(); ++iter)
							{
								LENGTH temp = TabPhi[iter->eid1]+TabPhi[iter->eid2];
								if(temp < newphie)
									newphie = temp;
							}
						}
					}

					bool eNotInQ=QQ.notIn(e);
					//assert(Q.index[e]!=Q.NO);//e should not be pop out before


					if(e < iENum0)
						newphie = min(eh.edges[e].length, newphie);
					if(bk == newphie)
						continue;
					if(eNotInQ)
					{
						TabPhi[e] = newphie;
						if(pi[c]<pi[y])
						{
							KEY_TYPE ktemp(c,y,TabPhi[e]);
							QQ.Push(ktemp, e);
							//if(!eh.existEdge(c, y, e))
						//assert(false);
						}
						else
						{
							KEY_TYPE ktemp(y,c,TabPhi[e]);
							QQ.Push(ktemp, e);
							//if(!eh.existEdge(c, y, e))
						//assert(false);
						}
					}
					else if(bk > newphie)
					{
						TabPhi[e] = newphie;
						//assert(Q.minheap[Q.index[e]].second == e);
						QQ.minheap[QQ.index[e]].first.key3 = newphie;
					}
				}
			}
		}
		static int timestamp=0;
		clock_t start = clock();
		updateDis(r, timestamp);
		dCaldistTime += clock()-start;
		++timestamp;
	}

	void update(int eid)
	{
		int u = eh.edges[eid].start;
		int v = eh.edges[eid].end;
		LENGTH ww = eh.edges[eid].length+2;
		update(eid, u, v, ww);
	}

	void updates(int count)
	{
		for(int ii=0; ii<count; ++ii)
		{
			int eid = rand() / (double)RAND_MAX*iENum0;
			if(eid>-1 && eid<iENum0)
				update(eid);
		}
	}

	void batchs(const vector<pair<int, LENGTH> > & vBatch)
	{
		KEY_TYPE::s_pi = &pi;
		vector<KEY_TYPE> vForSort;
		for(int ii=0; ii<vBatch.size(); ++ii)
		{
			int eid = vBatch[ii].first;
			int u = eh.edges[eid].start;
			int v = eh.edges[eid].end;
			if(pi[u]>pi[v])
				swap(u,v);

			KEY_TYPE temp(u,v,ii);//ii is not the phi value here,just for index
			vForSort.push_back(temp);
		}
		sort(vForSort.begin(), vForSort.end());

		MinHeap2 Q;
		
		Q.initial(iENum);
		

		for(int ii=0; ii < vForSort.size(); ++ii)
		{
			int jj = vForSort[ii].key3;
			int eid = vBatch[jj].first;
			LENGTH ww = vBatch[jj].second;
			LENGTH newphi=(LENGTH)INFINITY;
			if(ww < TabPhi[eid])
				newphi = ww;
			else
			{
				const set<TRIPLE> & triple = TabE[eid];
				LENGTH newphi = (int)INFINITY;
				for(set<TRIPLE>::iterator iter=triple.begin(); iter!=triple.end(); ++iter)
				{
					LENGTH temp = TabPhi[iter->eid1]+TabPhi[iter->eid2];
					if(temp < newphi)
						newphi = temp;
				}
				newphi=min(ww, newphi);
				if(newphi == TabPhi[eid])
				{
					eh.edges[eid].length = ww;
					continue;
				}
			}
			TabPhi[eid]=newphi;
			eh.edges[eid].length=ww;
			KEY_TYPE ktemp(vForSort[ii].key1,vForSort[ii].key2,newphi);
			Q.Push(ktemp,eid);
		}

		vector<int> vR;
		while(Q.count)
		{
			vector<pair<KEY_TYPE, int> > vQ;
			vQ.push_back(Q.minheap[0]);

			int eid;
			KEY_TYPE currentKey = Q.Pop(eid);
			TabPhi[eid]=currentKey.key3;
			int v1 = currentKey.key1;
			int r = v1;
		
			while(Q.count)
			{
				if(Q.minheap[0].first.key1 == v1)
				{
					vQ.push_back(Q.minheap[0]);
					int temp;
					KEY_TYPE ktemp=Q.Pop(temp);
					TabPhi[temp] = ktemp.key3;
				}
				else
					break;
			}
			int iSizevQ=vQ.size();
			TreeNode & tn=TG[v1];
			vector<int> & vNeighbor = tn.neighbors;
			bool bNoUpdate=true;
			for (int ii = 0; ii < vNeighbor.size(); ++ii)
			{
				int c = vNeighbor[ii];
				//int e1=-1;
				//eh.existEdge(v1, c, e1);
				//LENGTH phi1 = TabPhi[e1];

				for(int jj = 0; jj < iSizevQ; ++jj)
				{
					int y=vQ[jj].first.key2;
					if(c==y)
						continue;
					//int e2=-1;
					//eh.existEdge(v1, y, e2);
					int e=-1;
					eh.existEdge(c, y, e);

					//should consider whether e is in Q or not
					//if not, if phi remain, ignore, else update.
					//if yes, if phi increase, or remian, ignore; else update phi.
					LENGTH bk=TabPhi[e];
					bool eNotInQ=Q.notIn(e);
					//assert(Q.index[e]!=Q.NO);//e should not be pop out before

					const set<TRIPLE> & triple = TabE[e];
					LENGTH newphie=(LENGTH)INFINITY;
					for(set<TRIPLE>::iterator iter=triple.begin(); iter!=triple.end(); ++iter)
					{
						LENGTH temp = TabPhi[iter->eid1]+TabPhi[iter->eid2];
						if(temp < newphie)
							newphie = temp;
					}
					if(e < iENum0)
						newphie = min(eh.edges[e].length, newphie);
					if(bk == newphie)
						continue;
					if(eNotInQ)
					{
						TabPhi[e] = newphie;
						bNoUpdate = false;
						if(pi[c]<pi[y])
						{
							KEY_TYPE ktemp(c,y,TabPhi[e]);
							Q.Push(ktemp, e);
							//if(!eh.existEdge(c, y, e))
						//assert(false);
						}
						else
						{
							KEY_TYPE ktemp(y,c,TabPhi[e]);
							Q.Push(ktemp, e);
							//if(!eh.existEdge(c, y, e))
						//assert(false);
						}
					}
					else if(bk > newphie)
					{
						TabPhi[e] = newphie;
						assert(Q.minheap[Q.index[e]].second == e);
						Q.minheap[Q.index[e]].first.key3 = newphie;
						bNoUpdate = false;
					}
				}
			}
			if(bNoUpdate)
				vR.push_back(r);
		}
		//sort vR by pi, remove duplicate.
		vForSort.resize(0);
		for(int ii=0; ii<vR.size(); ++ii)
		{
			KEY_TYPE temp(vR[ii], 0, 0);//0,0 is not used at all
			vForSort.push_back(temp);
		}
		sort(vForSort.begin(), vForSort.end());
		int rSize = vR.size();
		static int timestamp = 1;
		clock_t start=clock();
		for(int ii=rSize-1; ii>-1; --ii)
		{
			int r = vForSort[ii].key1;
			updateDis(r, timestamp);//duplicate will be handled by timestamp, timestamp handle more than that.
		}
		dCaldistTime += clock()-start;
		++timestamp;
	}

	void batchs()
	{
		clock_t start=clock();
		int iSize = vUpdatedEdges.size();
		batchs(vUpdatedEdges);
		double dtime = (double)(clock() - start) / CLOCKS_PER_SEC;
		cout <<"batch update time="<<dtime<< endl;
		cout <<"avg time="<<dtime/iSize<<endl;
		cout <<"edges="<<iSize<<endl;
	}


	void getDescendants(int a, vector<int> &vD)
	{
		if(a<0)
			return;
		vD.push_back(a);
		for(int ii=0; ii < TG[a].children.size(); ++ii)
		{
			getDescendants(TG[a].children[ii], vD);
		}
	}

	vector<bool> vUpdated;
	vector<pair<int,int> > vTestPair;

	vector<pair<int, LENGTH> > vUpdatedEdges;
	vector<pair<int, int> > vUpdatedVertices;

	void writeUpdateEdgePrepareTest(ofstream & fU, int count = 1)
	{
		for(int no=0; no < count; ++no){
		int eid,u,v;
		while(true)
		{
			eid = rand() / (double)RAND_MAX*iENum0;
			if(eid>=iENum0 ||  vUpdated[eid])
				continue;

			u = eh.edges[eid].start;
			v = eh.edges[eid].end;
			if(pi[u]>pi[v])
				swap(u,v);

			//check if node u is a branch node
			const int iNum = TG[u].children.size() ;
			if(iNum>1)
				break;
		}
		
		const int h1=TG[u].height;
		const int h2=TG[v].height;
		//random generate query involves the update edge
		//find the desecndants of node u
		//get the first two children of u
		int c1 = TG[u].children[0];
		int c2 = TG[u].children[1];

		vector<int> vD1, vD2;
		getDescendants(c1, vD1);
		getDescendants(c2, vD2);

		vector<pair<int,int> > vPairs;
		vector<LENGTH> vDisOld;

		int s1 = vD1.size(), s2 = vD2.size();
		LENGTH d1,d2, dmin;
		for(int ii=0; ii < s1; ++ii)
		{
			int v1=vD1[ii];
			for(int jj = 0; jj < s2; ++jj)
			{
				int v2=vD2[jj];
				d1 = TG[v1].vToAncestor[h1]+TG[v2].vToAncestor[h1];
				d2 = TG[v1].vToAncestor[h2]+TG[v2].vToAncestor[h2];
				if(d1 != d2)
					continue;
				dmin=d1;
				bool bshort=true;
				for(int tt=0; tt<TG[u].pos.size(); ++tt)
				{
					LENGTH d = TG[v1].vToAncestor[TG[u].pos[tt]]+TG[v2].vToAncestor[TG[u].pos[tt]];
					if(d < dmin)
					{
						bshort = false;
						break;
					}
				}
				if(bshort)
				{
					//find a pair
					pair<int,int> temp(v1, v2);
					vPairs.push_back(temp);
					vDisOld.push_back(dmin);
				}

			}
		}
		if(vPairs.size()==0)
		{
			//checkCorre
			continue;
		}

		//do update	
		vUpdated[eid] = true;
		LENGTH ww = eh.edges[eid].length*2*(rand() / (double)RAND_MAX);
		//cout<<"\nedge "<<eid<<"=("<<u<<","<<v<<"),weight:"<<eh.edges[eid].length<<" to "<<ww<<endl;
		fU<<eid<<'\t'<<u<<'\t'<<v<<'\t'<<ww<<'\t'<<ww-eh.edges[eid].length<<endl;

		vTestPair.insert(vTestPair.end(),vPairs.begin(), vPairs.end());

		pair<int, int> tempv(u,v);
		vUpdatedVertices.push_back(tempv);
		pair<int, LENGTH> tempe(eid, ww);
		vUpdatedEdges.push_back(tempe);

		}
	}

	void writeTestPair(ofstream & fT)
	{
		int iSize = vTestPair.size();
		for(int ii=0; ii<iSize; ++ii)
		{
			int v1=vTestPair[ii].first;
			int v2=vTestPair[ii].second;
			LENGTH dis=(LENGTH)INFINITY;

			int vc=lca.LCAQuery(v1,v2);
			if(vc == v1)
			{
				dis = TG[v2].vToAncestor[TG[v1].height];
			}
			else if(vc == v2)
			{
				dis = TG[v1].vToAncestor[TG[v2].height];
			}
			else
			{
				int ps = TG[vc].pos.size();
				int * p2 = TG[vc].pos.data();
				LENGTH *dx = TG[v1].vToAncestor.data(), *dy = TG[v2].vToAncestor.data();
				for (int i = 0; i < ps; i++){
					LENGTH tmp = dx[p2[i]] + dy[p2[i]];
					if (dis > tmp)
						dis = tmp;
				}
			}
			fT<<ii<<'\t'<<v1<<'\t'<<v2<<'\t'<<dis<<endl;
		}
	}

	void checkTest()
	{
		string strTests("testPairs.txt");
		string strFail("failCase.txt");
		ifstream fT(strTests.c_str());
		ofstream fF(strFail.c_str());
		while(!fT.eof())
		{
			int v1,v2,id;
			LENGTH dis;
			fT>>id>>v1>>v2>>dis;
			LENGTH res=(LENGTH)INFINITY;
			int vc=lca.LCAQuery(v1,v2);
			if(vc == v1)
			{
				res = TG[v2].vToAncestor[TG[v1].height];
			}
			else if(vc == v2)
			{
				res = TG[v1].vToAncestor[TG[v2].height];
			}
			else
			{
				int ps = TG[vc].pos.size();
				int * p2 = TG[vc].pos.data();
				LENGTH *dx = TG[v1].vToAncestor.data(), *dy = TG[v2].vToAncestor.data();

				//LENGTH res = (int)INFINITY;
				for (int i = 0; i < ps; i++){
					LENGTH tmp = dx[p2[i]] + dy[p2[i]];
					if (res > tmp)
						res = tmp;
				}
				
			}
			if(res != dis)
				fF<<id<<'\t'<<v1<<'\t'<<v2<<'\t'<<dis<<'\t'<<res<<endl;
		}
		fT.close();
		fF.close();
	}

	void checkTestBatch()
	{
		string strUpdate("updateEdges.txt");
		string strTests("testPairs.txt");
		string strFail("failCaseBatch.txt");
		ifstream fU(strUpdate.c_str());
		ifstream fT(strTests.c_str());
		ofstream fF(strFail.c_str());

		vector<pair<int, LENGTH> > vBatch;

		while(!fU.eof())
		{
			pair<int, LENGTH> temp;
			int u,v;
			LENGTH delta;
			fU>>temp.first>>u>>v>>temp.second>>delta;
			vBatch.push_back(temp);
		}
		fU.close();

		batchs(vBatch);

		while(!fT.eof())
		{
			int v1,v2,id;
			LENGTH dis;
			fT>>id>>v1>>v2>>dis;
			LENGTH res=(LENGTH)INFINITY;
			int vc=lca.LCAQuery(v1,v2);
			if(vc == v1)
			{
				res = TG[v2].vToAncestor[TG[v1].height];
			}
			else if(vc == v2)
			{
				res = TG[v1].vToAncestor[TG[v2].height];
			}
			else
			{
				int ps = TG[vc].pos.size();
				int * p2 = TG[vc].pos.data();
				LENGTH *dx = TG[v1].vToAncestor.data(), *dy = TG[v2].vToAncestor.data();

				//LENGTH res = (int)INFINITY;
				for (int i = 0; i < ps; i++){
					LENGTH tmp = dx[p2[i]] + dy[p2[i]];
					if (res > tmp)
						res = tmp;
				}
				
			}
			if(res != dis)
				fF<<id<<'\t'<<v1<<'\t'<<v2<<'\t'<<dis<<'\t'<<res<<endl;
		}
		fT.close();
		fF.close();
	}

	void updateEdgeOBO()
	{
		QQ.initial(iENum);
		clock_t start=clock();
		int iSize = vUpdatedEdges.size();
		for(int ii = 0; ii < iSize; ++ii)
		{
			int eid = vUpdatedEdges[ii].first;
			LENGTH ww = vUpdatedEdges[ii].second;
			int u = vUpdatedVertices[ii].first;
			int v = vUpdatedVertices[ii].second;

			update(eid,u,v,ww);
		}
		double dtime = (double)(clock() - start) / CLOCKS_PER_SEC;
		cout <<"one by one update time="<<dtime<< endl;
		cout <<"avg time="<<dtime/iSize<<endl;
		cout <<"edges="<<iSize<<endl;
	}

	void updateRandOBO(const vector<double> & vRatio, int iCount)
	{
		QQ.initial(iENum);
		vUpdated.resize(iENum0, false);
		int num = 0;
		vector<int> vUpdatedEdges;
		vector<pair<int, int> > vUpdatedVertices;
		while(num < iCount)
		{
			int eid = rand()*1.0/RAND_MAX * iENum0;
			if(eid < iENum0 && vUpdated[eid] == false)
			{
				vUpdated[eid] = true;
				int u = eh.edges[eid].start;
				int v = eh.edges[eid].end;
				pair<int, int> tmp(u,v);
				vUpdatedVertices.push_back(tmp);
				vUpdatedEdges.push_back(eid);
				++num;
			}
		}

		cout << endl;
		for(int ii = 0; ii < vRatio.size(); ++ii)
		{
			//cout << ii <<':';
			vector<LENGTH> vUpdatedWeight(iCount);
			for(int jj = 0; jj < iCount; ++jj)
			{
				vUpdatedWeight[jj] = eh.edges[vUpdatedEdges[jj]].length*(1+vRatio[ii]);
				if(vUpdatedWeight[jj]<0)
					vUpdatedWeight[jj] = 0;
			}
			clock_t start=clock();
			for(int jj = 0; jj < iCount; ++jj)
				update(vUpdatedEdges[jj], vUpdatedVertices[jj].first, vUpdatedVertices[jj].second, vUpdatedWeight[jj]);
			double davgtime = (double)(clock() - start) /CLOCKS_PER_SEC/  iCount;
			cout << davgtime <<'\t';
		}
		cout << endl;
	}

	void testDec(int iCount = 1000)
	{
		int iGroup = 1;
		vector<double> vRatio(iGroup);
		for(int ii = 0; ii < iGroup; ++ii)
			vRatio[ii] = -rand()*1.0/RAND_MAX;
		dCaldistTime = 0.0;
		updateRandOBO(vRatio, iCount);
		cout << dCaldistTime * 1.0/CLOCKS_PER_SEC / iCount<<endl;
	}
	void testInc(int iCount = 1000)
	{
		int iGroup = 1;
		vector<double> vRatio(iGroup);
		for(int ii = 0; ii < iGroup; ++ii)
			vRatio[ii] = rand()*1.0/RAND_MAX;
		dCaldistTime = 0.0;
		updateRandOBO(vRatio, iCount);
		cout << dCaldistTime * 1.0/CLOCKS_PER_SEC/iCount<<endl;
	}

	void updateRandBatch(const vector<double> & vCountPercentage, double ratio)
	{
		vector<int> vIDs(iENum0);
		for(int ii = 0; ii < iENum0; ++ii)
			vIDs[ii] = ii;
		random_shuffle(vIDs.begin(), vIDs.end());

		QQ.initial(iENum);
		vector<pair<int, LENGTH> > vBatch;
		int iGroup = vCountPercentage.size();
		int iStart = 0;
		for(int ii = 0; ii < iGroup; ++ii)
		{
			int iCount = iENum0 * vCountPercentage[ii];
			if(iCount <= 0)
				iCount = 1;
			if(iCount > iENum0)
				iCount = iENum0;
			for(int jj = iStart; jj < iCount; ++jj)
			{
				int eid = vIDs[jj];
				LENGTH ww = eh.edges[eid].length * ratio;
				pair<int, LENGTH> tmp(eid, ww);
				vBatch.push_back(tmp);
			}
			iStart = iCount;
			clock_t start=clock();
			batchs(vBatch);
			double davgtime = (double)(clock() - start) / iCount;
			cout << davgtime << '\t';
		}
		cout << endl;
	}

	void testDecBatch()
	{
		int iGroup = 1;
		vector<double> vCountPercentage(iGroup);
		for(int ii = 0; ii < iGroup; ++ii)
			vCountPercentage[ii] = 1000.0/iENum0;//pow(10.0, ii-iGroup);
		dCaldistTime = 0.0;
		updateRandBatch(vCountPercentage, 0.9);
		cout << dCaldistTime * 1.0/vCountPercentage[0]/iENum0<<endl;
	}

	void testIncBatch()
	{
		int iGroup = 1;
		vector<double> vCountPercentage(iGroup);
		for(int ii = 0; ii < iGroup; ++ii)
			vCountPercentage[ii] = 1000.0/iENum0;//pow(10.0, ii-iGroup);
		dCaldistTime = 0.0;
		updateRandBatch(vCountPercentage, 1.1);
		cout << dCaldistTime * 1.0/vCountPercentage[0]/iENum0<<endl;
	}

	void testAll()
	{
		calSize();
		
		//testInc();
		
		testIncBatch();
		testDecBatch();
		//testDec();
		
	}

	void calSize()
	{
		double dSize1 = TabPhi.size()*sizeof(LENGTH);
		int iSize = TabE.size();
		double dSize2 = 0;
		for(int ii = 0; ii < iSize; ++ii)
		{
			const set<TRIPLE> & triple = TabE[ii];
			dSize2 += triple.size();
		}
		dSize2 *= sizeof(TRIPLE);
		cout<< dSize1/1024/1024 << "\tMB"<<endl;
		cout<< dSize2/1024/1024 << "\tMB"<<endl;
	}

};
double Node::s_outputThredhold = 0.5;


int main1(int argc, char* argv[])
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
	clock_t start = clock();

	ReadFile rf(filename);

	rf.readRoadNet(filename);

	char indexName[100];
	sprintf(indexName, "%s_%d-%d-%d.pindex", dataset.c_str(), method, bRoLCA, bheight_control);

	//rf.edgeList.component(47437, rf.vertexList);

	Tree_Decomposition td(rf,dataset);
	td.g_bRoLCA = bRoLCA;
	td.run(method);	
	td.calDisAndPos(td.root);
	td.lca.makeRMQ(td.root);
	cout << endl << "build time=" << (double)(clock() - start) / CLOCKS_PER_SEC  << endl;

	td.vUpdated.resize(td.iENum0, false);
	//td.update(256511,186605,186613,1042);
	string strUpdate(dataset+"updateEdges.txt");
	ofstream fU(strUpdate.c_str());//,std::ofstream::app);
	td.writeUpdateEdgePrepareTest(fU, method);
	fU.close();

	string strTests("testPairs.txt");
	ofstream fT(strTests.c_str());
	td.writeTestPair(fT);
	fT.close();

	if(bRoLCA)
		td.batchs();
	else
		td.updateEdgeOBO();
	//td.batchs();

	td.write(dataset+".index");
	td.writeTable(dataset+".table");
	return 0;
}

int main2(int argc, char* argv[])
{
	string dataPath("../data/");
	string dataset("NY");

	string strUpdate("updateEdges.txt");
	string strTests("testPairs.txt");

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
	clock_t start = clock();

	ReadFile rf(filename);

	rf.readRoadNet(filename);
	rf.updateEdges(strUpdate.c_str());

	char indexName[100];
	sprintf(indexName, "%s_%d-%d-%d.pindex", dataset.c_str(), method, bRoLCA, bheight_control);


	//rf.edgeList.component(47437, rf.vertexList);

	Tree_Decomposition td(rf);


	td.g_bRoLCA = bRoLCA;
	td.run(method);
	cout << (double)(clock() - start) / CLOCKS_PER_SEC  << '\t';
	td.resetHeight(td.root);
	td.calDisAndPos(td.root);

	
	td.avgHeight();
	td.avgWidth();
	td.lca.makeRMQ(td.root);

	td.checkTest();

	return 0;
}

int main3(int argc, char* argv[])
{
	string dataPath("../data/");
	string dataset("NY");

	string strUpdate("updateEdges.txt");
	string strTests("testPairs.txt");

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
	clock_t start = clock();

	ReadFile rf(filename);

	rf.readRoadNet(filename);
	rf.updateEdges(strUpdate.c_str());

	char indexName[100];
	sprintf(indexName, "%s_%d-%d-%d.pindex", dataset.c_str(), method, bRoLCA, bheight_control);


	//rf.edgeList.component(47437, rf.vertexList);

	Tree_Decomposition td(rf);


	td.g_bRoLCA = bRoLCA;
	td.run(method);
	cout << (double)(clock() - start) / CLOCKS_PER_SEC  << '\t';
	td.resetHeight(td.root);
	td.calDisAndPos(td.root);

	
	td.avgHeight();
	td.avgWidth();
	td.lca.makeRMQ(td.root);

	td.checkTestBatch();

	return 0;
}

int main4(int argc, char* argv[])
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
	string filename(dataPath + dataname);
	clock_t start = clock();
	ReadFile rf(filename);
	rf.readRoadNet(filename);

	Tree_Decomposition td(rf);

	td.g_bRoLCA = bRoLCA;
	td.run(method);
	cout << (double)(clock() - start) / CLOCKS_PER_SEC  << '\t';
	td.resetHeight(td.root);
	td.calDisAndPos(td.root);	
	td.avgHeight();
	td.avgWidth();
	td.lca.makeRMQ(td.root);
 
	td.testAll();

	return 0;
}

int main(int argc, char* argv[])
{
	char what=argv[5][0];
	if(what=='i')
		main1(argc, argv);
	else if(what=='t')
		main2(argc, argv);
	else if(what=='b')
		main3(argc, argv);
	else
		main4(argc, argv);

	return 0;
}