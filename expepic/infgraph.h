#include "iheap.h"
#include <queue>	//priority_queue
#include <utility>  // pair
#include <numeric>

class InfGraph: public Graph
{
private:
    vector<bool> visit;
    vector<int> visit_mark;
	vector<int> curNode;	
	vector<int> curNodeIdx;
public:
	unsigned int NumcurNode;

    vector<vector<int>> hyperG;
    vector<vector<int>> hyperGT;
	
	vector<vector<int>> hyperG_2;
	vector<vector<int>> hyperGT_2;

	vector<bool>active_set;
	vector<vector<int>>PO;
		
	sfmt_t sfmtSeed;
	vector<int> seedSet;	

    InfGraph(string folder, string graph_file): Graph(folder, graph_file)
    {
        srand(time(NULL));
        sfmt_init_gen_rand(&sfmtSeed , rand());		

        visit = vector<bool> (n);
        visit_mark = vector<int> (n);
		
		curNode = vector<int>(n);
		curNodeIdx = vector<int>(n);
		iota(curNode.begin(), curNode.end(), 0);
		iota(curNodeIdx.begin(), curNodeIdx.end(), 0);
		NumcurNode = n;

		active_set = vector<bool>(n, false);				

		hyperG.resize(n, vector<int>());		

		hyperG_2.resize(n, vector<int>());

		PO.resize(n, vector<int>());
    }

    void init_hyper_graph(){

		for (auto& hyper : hyperG)hyper.clear();
		for (auto& hyperT : hyperGT)vector<int>().swap(hyperT);
		hyperGT.clear();

		for (auto& hyper : hyperG_2)hyper.clear();
		for (auto& hyperT : hyperGT_2)vector<int>().swap(hyperT);
		hyperGT_2.clear();
		
		seedSet.clear(); 		
		fill(active_set.begin(), active_set.end(), false);
		
		iota(curNode.begin(), curNode.end(), 0);
		iota(curNodeIdx.begin(), curNodeIdx.end(), 0);
		NumcurNode = n;
    }

	//initialize PO each time it loads a possible world.
	void load_possible_world(string index, const Argument & arg)
	{		
		//initialization
		for (unsigned int i = 0; i < PO.size(); ++i)PO[i].clear();

		int index1 = 0;
		while (arg.dataset[index1] != '/')index1++;
		int index2 = ++index1;
		while (arg.dataset[index2] != '/')index2++;

		string file_name = arg.dataset + arg.dataset.substr(index1, index2 - index1) + "_" + index;

		if (influModel == LT)file_name += "_lt";		

		size_t length;
		auto ptr = map_file(file_name.c_str(), length);
		auto f = ptr;
		auto tep = f;
		char buffer[256];

		auto l = f + length;

		unsigned int t1, t2;

		while (f && f != l)
		{
			if ((f = static_cast<char*>(memchr(f, '\n', l - f))))
			{
				memset(buffer, 0, sizeof(buffer));
				memcpy(buffer, tep, f - tep);
				f++;
				tep = f;
				stringstream s;
				s << buffer;
				s >> t1 >> t2;
				ASSERT(t1 < n);
				ASSERT(t2 < n);

				PO[t1].push_back(t2);												
			}
		}
		munmap(ptr, length);
	}

	char* map_file(const char* fname, size_t& length)
	{
		int fd = open(fname, O_RDONLY);
		if (fd == -1)
			handle_error("open");

		// obtain file size
		struct stat sb;
		if (fstat(fd, &sb) == -1)
			handle_error("fstat");

		length = sb.st_size;

		char* addr = static_cast<char*>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));
		if (addr == MAP_FAILED)
			handle_error("mmap");

		// TODO close fd at some point in time, call munmap(...)
		close(fd);
		return addr;
	}

    void build_hyper_graph_r(unsigned int R)
    {
        if( R > INT_MAX ){
            cout<<"Error:R too large"<<endl;
            exit(1);
        }    
		//every time starts from scratch
		auto counter = hyperGT.size();
		while (counter < R)BuildHypergraphNode(counter++);		
    }	

    #include "discrete_rrset.h"

	double build_seedset(int targetSize, vector<int>&batch_set)
	{
		
		vector<int> coverage(n, 0);
		size_t maxDeg = 0;
		for (auto i = n; i--;)
		{
			//const auto deg = node_influence_1[i];
			const auto deg = hyperG[i].size();
			coverage[i] = deg;
			if (deg > maxDeg) maxDeg = deg;
		}
		vector<vector<int>> degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree
		for (auto i = n; i--;)
		{
			//if (coverage[i] == 0) continue;
			degMap[coverage[i]].push_back(i);
		}
		vector<int> sortedNode(n); // sortedNode: record the sorted nodes in ascending order of degree
		vector<int> nodePosition(n); // nodePosition: record the position of each node in the sortedNode
		vector<int> degreePosition(maxDeg + 2); // degreePosition: the start position of each degree in sortedNode
		uint32_t idxSort = 0;
		size_t idxDegree = 0;
		for (auto& nodes : degMap)
		{
			degreePosition[idxDegree + 1] = degreePosition[idxDegree] + (uint32_t)nodes.size();
			idxDegree++;
			for (auto& node : nodes)
			{
				nodePosition[node] = idxSort;
				sortedNode[idxSort++] = node;
			}
		}
		// check if an edge is removed
		std::vector<bool> edgeMark(hyperGT.size(), false);
		// record the total of top-k marginal gains
		int sumTopk = 0;
		for (auto deg = maxDeg + 1; deg--;)
		{
			if (degreePosition[deg] <= (int)(n - targetSize))
			{
				sumTopk += deg * (degreePosition[deg + 1] - (int)(n - targetSize));
				break;
			}
			sumTopk += deg * (degreePosition[deg + 1] - degreePosition[deg]);
		}
		double boundMin = 1.0 * sumTopk;		
		double sumInf = 0;

		batch_set.clear();
		
		/*
		* sortedNode: position -> node
		* nodePosition: node -> position
		* degreePosition: degree -> position (start position of this degree)
		* coverage: node -> degree
		* e.g., swap the position of a node with the start position of its degree
		* swap(sortedNode[nodePosition[node]], sortedNode[degreePosition[coverage[node]]])
		*/
		for (auto k = targetSize; k--;)
		{
			const auto seed = sortedNode.back();
			sortedNode.pop_back();
			int newNumV = sortedNode.size();
			sumTopk += coverage[sortedNode[newNumV - targetSize]] - coverage[seed];			
			sumInf += coverage[seed];
			batch_set.push_back(seed);
			coverage[seed] = 0;
			for (auto edgeIdx : hyperG[seed])
			{
				if (hyperGT[edgeIdx].size() == 0 || edgeMark[edgeIdx]) continue;
				edgeMark[edgeIdx] = true;
				for (auto nodeIdx : hyperGT[edgeIdx])
				{
					if (coverage[nodeIdx] == 0) continue; // This node is seed, skip
					const auto currPos = nodePosition[nodeIdx]; // The current position
					const auto currDeg = coverage[nodeIdx]; // The current degree
					const auto startPos = degreePosition[currDeg]; // The start position of this degree
					const auto startNode = sortedNode[startPos]; // The node with the start position
					// Swap this node to the start position with the same degree, and update their positions in nodePosition
					std::swap(sortedNode[currPos], sortedNode[startPos]);
					nodePosition[nodeIdx] = startPos;
					nodePosition[startNode] = currPos;
					// Increase the start position of this degree by 1, and decrease the degree of this node by 1
					degreePosition[currDeg]++;
					coverage[nodeIdx]--;
					// If the start position of this degree is in top-k, reduce topk by 1
					if (startPos >= newNumV - targetSize) sumTopk--;
				}
			}
			const auto boundLast = 1.0 * (sumInf + sumTopk);
			if (boundMin > boundLast) boundMin = boundLast;			
		}
		//cout << batch_set.size() << ": boundMin: " << boundMin << endl;// " boundLast: " << boundLast << endl;// << " sumTopk: " << sumTopk << " sum: " << sumInf + sumTopk << endl;
		return 1.0*boundMin;
	}
 
	double Coverage(unsigned int sample, vector<int>& batch_set)
	{
		auto idx = hyperGT_2.size();
		while (idx < sample)BuildHypergraph(idx++);
		vector<bool>overlap(idx, false);
		for (auto seed : batch_set)
			for (auto rrIdx : hyperG_2[seed])overlap[rrIdx] = true; //rrIdx is the index of the rr set that contains seedIdx.					
		
		return count(overlap.begin(), overlap.end(), true);
	}
};