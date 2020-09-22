#ifndef READFILE_H
#define READFILE_H

#include "Graph.h"
class ReadFile
{
public:
	bool m_bRecordServerOnEdge;
	string dataFile;
	string coordinateFile;
	string selectEdgeFile;
	ReadFile(string datafile);

	void readRoadNet(string filename);
	void pruneEdges();
	void readData();

	void groupForDegreeTwo(int vid, vector<set<int> > & vsGroup, vector<int> & vGroup, int direction);

	void readCoordinate(string filename);

	void writeMap();

	LENGTH singleSource(int source, int destination);

	int m_nVertex;
	int m_nEdge;

	EdgeList edgeList;
	VertexList vertexList;
private:
	void readHead(ifstream & fin);
};

#endif
