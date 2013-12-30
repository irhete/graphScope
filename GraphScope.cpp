#include <iostream>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <map>
#include <queue>
#include <bitset>
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

struct GraphSegment {
	string adjacency_matrix;

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
		map<int, int> &nodes, map<int, int> &renumbered_nodes) {
	string line;
	ifstream myfile(filename.c_str());
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
				renumbered_nodes[node_index] = v1;
				node_index++;
			}
			iter = nodes.find(v2);
			if (iter == nodes.end()) {
				nodes[v2] = node_index;
				renumbered_nodes[node_index] = v2;
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
	map<int, int> renumbered_nodes;
	const string filename = "small_graph.txt";

	read_graph(filename, edges, nodes, renumbered_nodes);
	cout << nodes.size() << endl;
	int vertex_count = 4;
	const int bitset_size = 16;
	vector<bitset<bitset_size> > segments;

	bitset<bitset_size> segment;
	while (!edges.empty()) {
		Edge e = edges.top();
		while (!edges.empty() && edges.top().timestamp == e.timestamp) {
			e = edges.top();
			edges.pop();
			cout << e.v1 << ", " << e.v2 << ", " << e.timestamp << endl;
			segment.set(vertex_count * e.v1 + e.v2);
			segment.set(vertex_count * e.v2 + e.v1);
		}
		cout << segment << endl;
		segments.push_back(segment);
	}

	return 0;

}
