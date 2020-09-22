#include "Readfile.h"
#include <fstream>

#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <string>
#include "Heap.h"

ReadFile::ReadFile(string datafile)
	:dataFile(datafile)
{
}

void ReadFile::readHead(ifstream & fin)
{
	int iLine = 0;
	string aline;
	while (iLine < 4)
	{
		getline(fin, aline);
		++iLine;
		//cout << aline << endl;
	}

	string temp;
	fin >> temp; fin >> temp;
	fin >> m_nVertex >> m_nEdge;
#ifdef SHOW_DETAIL
	cout << m_nVertex << '\t' << m_nEdge;
#endif

	while (iLine < 7)
	{
		getline(fin, aline);
		++iLine;
		//cout << aline << endl;
	}
}

void ReadFile::readRoadNet(string filename)
{
	//cout<<"Reading road network..."<<endl;
	ifstream fin(filename.c_str());
	if (!fin)
	{
		cout << "Readfile failed" << endl;
		return;
	}
	else
		readHead(fin);

	vertexList.resize(m_nVertex);

	Edge newEdge;
	Vertex newVertex;

	for (int ii = 0; ii < m_nVertex; ii++)
		vertexList[ii].id = ii;

	int edgeID=0;
	int iVertex = 0;
	char skipme;
	int iStart, iEnd;
	LENGTH length;
	while (!fin.eof() && edgeID < m_nEdge)
	{
		fin >> skipme >> iStart >> iEnd >> length;

		iStart = iStart - 1;
		iEnd = iEnd - 1;
		
		newEdge.start = iStart;
		newEdge.end = iEnd;
		newEdge.length = length;
		newEdge.id = edgeID;

		bool bInsert = edgeList.insertAnEdge(newEdge);

		if (bInsert)
		{
			vertexList[iStart].insertANeighborEdge(edgeID);
			//vertexList[iStart].neighbor_list.push_back(iEnd);

			vertexList[iEnd].insertANeighborEdge(edgeID);
			//vertexList[iEnd].neighbor_list.push_back(iStart);

			++edgeID;
		}
	}
	fin.close();

	pruneEdges();

	#ifndef SHOW_DETAIL
	//cout << iStart << '\t' << iEnd;
	//cout<<"EdegNum: "<<edgeID<<'\t'<<"VertexNum: "<<m_nVertex<<'\t';
	cout<<edgeID<<'\t'<<m_nVertex<<'\t';

	int degree = 0;
	int maxDegree = 0, avgDegree = 0, countDegree2 = 0, countDegree3 = 0;
	for(int ii = 0; ii < m_nVertex; ++ii)
	{
		degree = vertexList[ii].edge_list.size();
		if(degree > maxDegree)
			maxDegree = degree;

		avgDegree += degree;
		if(degree == 2)
			++countDegree2;
				if(degree == 3)
			++countDegree3;
	}
	//cout << "maxDegree = "<<maxDegree <<"  avgDegree = "<<avgDegree*1.0/m_nVertex<<"  PercentageDegree2 = " <<countDegree2*1.0/m_nVertex <<endl;
	cout <<countDegree2*1.0/m_nVertex <<'\t'<<countDegree3*1.0/m_nVertex <<'\t';
	/*
	vector<int> vGroup(m_nVertex, -1);

	vector<set<int> > vsGroup;
	for(int ii=0; ii<m_nVertex; ++ii)
	{
		degree = vertexList[ii].edge_list.size();
		if(degree == 2 && vGroup[ii] == -1)
		{
			groupForDegreeTwo(ii, vsGroup, vGroup, -1);
		}
	}

	for(int ii = 0; ii < vsGroup.size(); ++ii)
	{
		if(vsGroup[ii].size() > 3)
			cout << ii <<'\t';
	}
	*/
	#endif
}

void ReadFile::groupForDegreeTwo(int vid, vector<set<int> > & vsGroup, vector<int> & vGroup, int direction)
{
	//cout << "vid0=" << vid <<'\t' << vGroup[vid] << '\t' << direction << endl;
	if(vertexList[vid].edge_list.size() != 2)
		return;
	int group = -1;
	if(direction == -1)
	{
		//should make a new group
		group = vsGroup.size();
		set<int> sGroup;

		sGroup.insert(vid);
		vGroup[vid] = group;

		vsGroup.push_back(sGroup);
	}
	else
	{
		group = vGroup[direction];
		vsGroup[group].insert(vid);
		vGroup[vid] = group;
	}

	for(set<int>::iterator iter = vertexList[vid].edge_list.begin() ; iter != vertexList[vid].edge_list.end(); ++iter)
	{
		int start = edgeList.edges[*iter].adjVertex(vid);
		if(start != direction)
			groupForDegreeTwo(start, vsGroup, vGroup, vid);
	}
}


void ReadFile::readCoordinate(string filename)
{
	cout<<"Reading coordinate..."<<endl;
	ifstream fin(filename.c_str());
	if (!fin)
	{
		cout << "Readfile failed" << endl;
		return;
	}
	else
	{
		int iLine = 0;
		string aline;
		while (iLine < 7)
		{
			getline(fin, aline);
			++iLine;
			//cout << aline << endl;
		}
	}

	int iCount = 0;
	char skipme;
	int vid = 0;
	LENGTH x,y;
	while (!fin.eof())
	{

		fin >> skipme >> vid >> x >> y;
		++iCount;
		if(iCount < 10)
		{
			cout << skipme << ' ' << vid << ' ';
			cout << x << '\t' << y << endl;
		}
		vertexList[vid-1].coordinate[0] = x;
		vertexList[vid-1].coordinate[1] = y;
	}

	fin.close();
}

void ReadFile::updateEdges(string filename)
{
	ifstream fin(filename.c_str());
	if(!fin)
	{
		cout<<"Read update edges failed."<<endl;
		return;
	}
	int eid,u,v;
	LENGTH ww,dw;
	while(!fin.eof())
	{
		fin>>eid>>u>>v>>ww>>dw;
		edgeList.edges[eid].length=ww;
	}
	fin.close();
}

void ReadFile::pruneEdges()
{
	//for each vertex, check whether their is an edge between two of its neighbors.
	return;

	int iNum = vertexList.size();

	vector<int> vStart, vEnd, vEdge;
	for(int ii = 0; ii < iNum; ++ii)
	{
		set<int> & sNeighbor = vertexList[ii].edge_list;

		for(set<int>::iterator jj = sNeighbor.begin(); jj != sNeighbor.end(); ++jj)
		{
			int eid1 = *jj;
			int vid1 = edgeList.edges[eid1].adjVertex(ii);

			LENGTH l1 = edgeList.edges[eid1].length;

			for(set<int>::iterator kk = jj; kk != sNeighbor.end(); ++kk)
			{
				if (kk == jj)
					continue;
				int eid2 = *kk;
				int vid2 = edgeList.edges[eid2].adjVertex(ii);
				
				int eid3;
				if(edgeList.existEdge(vid1, vid2, eid3))
				{
					LENGTH l2 = edgeList.edges[eid2].length;
					LENGTH l3 = edgeList.edges[eid3].length;

					if(l1 + l2 <= l3)
					{
						//remove eid3
						assert(vertexList[vid1].edge_list.find(eid3) != vertexList[vid1].edge_list.end());
						assert(vertexList[vid2].edge_list.find(eid3) != vertexList[vid2].edge_list.end());
						if(find(vEdge.begin(), vEdge.end(), eid3) != vEdge.end())
							continue;
						vStart.push_back(vid1);
						vEnd.push_back(vid2);
						vEdge.push_back(eid3);
						//cout<<vid1<<'\t'<<vid2<<'\t'<<eid3<<endl;
					}
				}
			}
		}
	}

	int iCount = vEdge.size();
#ifndef SHOW_DETAIL
	cout<< "Num of Prune Edges = " << iCount <<endl;
#endif
	for(int ii = 0; ii < iCount; ++ii)
	{
		//cout<<vStart[ii]<<'\t'<<vEnd[ii]<<'\t'<<vEdge[ii]<<endl;
		vertexList[vStart[ii]].deleteANeighborEdge(vEdge[ii]);
		vertexList[vEnd[ii]].deleteANeighborEdge(vEdge[ii]);

		pair<int, int> key1(vStart[ii], vEnd[ii]);
		pair<int, int> key2(vEnd[ii], vStart[ii]);
		edgeList.mapVToE.erase(key1);
		edgeList.mapVToE.erase(key2);
		edgeList.edges[vEdge[ii]].mark = -1;
	}
}

void ReadFile::readData()
{
	srand(int(time(0)));
	//cout<<"Step1: Read data..."<<endl;
	return;
}

void ReadFile::writeMap()
{
	char filename[100];
	//sprintf(filename, "C%dS%d", m_nCustomer,m_nServer);
	ofstream output(filename);
	int i;

	i = 0;
	while (1)
	{
		//output << edgeList[i].start << '\t' << edgeList[i].end << '\t' << edgeList[i].length << '\t';

		++i;
		if (i<m_nEdge)
			output << endl;
		else
			break;
	}
}