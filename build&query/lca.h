#include "Graph.h"
#include <xmmintrin.h>
class LCA{
public:
	LCA():pTG(NULL), iNodes(0){}
	~LCA()
	{
		if(LOG2)
			delete [] LOG2;
		if(LOGD)
			delete [] LOGD;
		if(toRMQ)
			delete [] toRMQ;
	}

	void initial(vector<TreeNode> & TG, int _iNodes)
	{
		iNodes = _iNodes;
		pTG = & TG;
		LOG2 = (int*)malloc(sizeof(int) * (iNodes * 2 + 10));
		LOGD = (int*)malloc(sizeof(int) * (iNodes * 2 + 10));
		int k = 0, j = 1;
		for (int i = 0; i < iNodes * 2 + 10; i++)
		{
			if (i > j * 2)
			{
				j *= 2;
				k++;
			}
			LOG2[i] = k;
			LOGD[i] = j;
		}
		//_mm_prefetch(LOG2, _MM_HINT_T0);
		//_mm_prefetch(LOGD, _MM_HINT_T0);
	}

	vector<TreeNode> * pTG;
	int iNodes;

	int * toRMQ;
	vector<int> EulerSeq;
	vector< vector<int> > RMQIndex;
	int * LOGD;
	int * LOG2;

	void makeRMQDFS(int p, int height){

		toRMQ[p] = EulerSeq.size();
		EulerSeq.push_back(p);
		for (int i = 0; i < (* pTG)[p].children.size(); i++){
			makeRMQDFS((* pTG)[p].children[i], height + 1);
			EulerSeq.push_back(p);
		}
	}

	void makeRMQ(int root){
		EulerSeq.clear();
		toRMQ = (int *) malloc(sizeof(int)*iNodes);
		//toRMQ.resize(iNodes);

		makeRMQDFS(root, 1);

		RMQIndex.clear();
		RMQIndex.push_back(EulerSeq);
		int m = EulerSeq.size();
		for (int i = 2, k = 1; i < m; i = i * 2, k++){
			vector<int> tmp;
			tmp.clear();
			tmp.resize(EulerSeq.size());
			for (int j = 0; j < m - i; j++){
				int x = RMQIndex[k - 1][j], y = RMQIndex[k - 1][j + i / 2];
				if ((* pTG)[x].height < (* pTG)[y].height)
					tmp[j] = x;
				else tmp[j] = y;
			}
			RMQIndex.push_back(tmp);
		}
	}

	int LCAQuery(int _p, int _q){

		int p = toRMQ[_p];
		int q = toRMQ[_q];
		if (p > q){
			int x = p;
			p = q;
			q = x;
		}
		int len = q - p + 1;
		int i = 1, k = 0;
		while (i * 2 < len){
			i *= 2;
			k++;
		}
		q = q - i + 1;
		if ((* pTG)[RMQIndex[k][p]].height < (* pTG)[RMQIndex[k][q]].height)
			return RMQIndex[k][p];
		else return RMQIndex[k][q]; 
	}

	inline	int LCAQuery2(int _p, int _q){
		int p = toRMQ[_p], q = toRMQ[_q];
		
		if (p > q){
			int x = p;
			p = q;
			q = x;
		}
		const int len = q - p + 1;
		
		const int i = LOGD[len], k = LOG2[len];
		
		q = q - i + 1;
		q=RMQIndex[k][q];
		p=RMQIndex[k][p];

		if ((* pTG)[p].height < (* pTG)[q].height)
			return p;
		else return q; 
	}
};
