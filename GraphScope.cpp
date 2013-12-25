#include <iostream>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <map>
#include <queue>
using namespace std;

vector<string> &split(const string &s, char delim, vector<string> &elems) {
	stringstream ss(s);
	string item;
	while (getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

vector<string> split(const string &s, char delim) {
	vector<string> elems;
	split(s, delim, elems);
	return elems;
}

struct Edge {
	int v1, v2, timestamp;
	Edge(int v1, int v2, int timestamp) :
			v1(v1), v2(v2), timestamp(timestamp) {
	}
};

class CompareEdges {
public:
	bool operator()(Edge& e1, Edge& e2) {
		if (e1.timestamp > e2.timestamp)
			return true;
		return false;
	}
};

void read_graph(string filename,
		priority_queue<Edge, vector<Edge>, CompareEdges> &edges,
		map<int, int> &nodes) {
	string line;
	ifstream myfile(filename);
	if (myfile.is_open()) {
		int node_index = 0;
		while (getline(myfile, line)) {
			vector<string> splitRow = split(line, ';');
			int v1 = atoi(splitRow[0].c_str());
			int v2 = atoi(splitRow[1].c_str());
			int timestamp = atoi(splitRow[3].c_str());
			map<int, int>::iterator iter = nodes.find(v1);
			if (iter == nodes.end()) {
				nodes[v1] = node_index;
				node_index++;
			}
			iter = nodes.find(v2);
			if (iter == nodes.end()) {
				nodes[v2] = node_index;
				node_index++;
			}
			edges.push(Edge(nodes[v1], nodes[v2], timestamp));
		}
	}
	myfile.close();
}

int main() {
	priority_queue<Edge, vector<Edge>, CompareEdges> edges;
	map<int, int> nodes;

	read_graph("/wrk/stacc-sna/data_extract/sample_network_65.txt", edges,
			nodes);

	int i = 0;
	cout << edges.top().timestamp << endl;
	while (!edges.empty()) {
		i++;
		edges.pop();
	}

	cout << "Edges: " << i << endl;
	cout << "Nodes: " << nodes.size() << endl;

	return 0;

}
