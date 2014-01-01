#include <iostream>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <map>
#include <queue>
#include <bitset>
#include <cmath>
using namespace std;

const int VERTEX_COUNT = 7;

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

priority_queue<Edge, vector<Edge>, CompareEdges> read_graph(string filename,
		map<int, int> &nodes, map<int, int> &renumbered_nodes) {
	priority_queue<Edge, vector<Edge>, CompareEdges> edges;
	string line;
	ifstream myfile(filename.c_str());
	if (myfile.is_open()) {
		int node_index = 0;
		while (getline(myfile, line)) {
			vector<string> splitRow = split(line, ';');
			int v1 = atoi(splitRow[0].c_str());
			int v2 = atoi(splitRow[1].c_str());
			int timestamp = atoi(splitRow[2].c_str());
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
	return edges;
}

double node_cost(int* partitioning, vector<int> partition_members,
		vector<int>* adjacency, int node, int partition) {
	const int k = 2; // partitions count

	int possible_edges_from_node = VERTEX_COUNT - 1;
	int partition_size = partition_members.size();
	int possible_edges_from_partition = partition_size * VERTEX_COUNT - 1;
	int node_neighbors_count = adjacency[node].size();
	int partition_neighbors_count = 0;

	double node_neighbors_in_partition[k] = { 0 };
	for (vector<int>::iterator it = adjacency[node].begin();
			it != adjacency[node].end(); ++it) {
		int neighbor = *it;
		node_neighbors_in_partition[partitioning[neighbor]]++;
	}

	double partition_neighbors_in_partition[k] = { 0 };
	for (vector<int>::iterator it = partition_members.begin();
			it != partition_members.end(); it++) {
		for (vector<int>::iterator it2 = adjacency[*it].begin();
				it2 != adjacency[*it].end(); it2++) {
			int neighbor = *it2;
			partition_neighbors_in_partition[partitioning[neighbor]]++;
			partition_neighbors_count++;
		}
	}

	double entropy = 0;
	for (int i = 0; i < k; i++) {
		if (node_neighbors_in_partition[i] > 0) {
			entropy += (node_neighbors_in_partition[i]
					/ possible_edges_from_node)
					* (log(partition_neighbors_in_partition[i]
									/ possible_edges_from_partition) / log(2));
		}
	}
	if (possible_edges_from_partition != partition_neighbors_count) {
	entropy += ((double) (possible_edges_from_node - node_neighbors_count)
			/ possible_edges_from_node)
			* (log( (double)
					(possible_edges_from_partition - partition_neighbors_count)
							/ possible_edges_from_partition) / log(2));
	}
	return entropy;
}

vector<bitset<VERTEX_COUNT * VERTEX_COUNT> > create_bitset(
		priority_queue<Edge, vector<Edge>, CompareEdges> edges) {
	const int bitset_size = VERTEX_COUNT * VERTEX_COUNT;
	vector<bitset<bitset_size> > segments;
	bitset<bitset_size> segment;
	while (!edges.empty()) {
		Edge e = edges.top();
		while (!edges.empty() && edges.top().timestamp == e.timestamp) {
			e = edges.top();
			edges.pop();
			cout << e.v1 << ", " << e.v2 << ", " << e.timestamp << endl;
			segment.set(VERTEX_COUNT * e.v1 + e.v2);
			segment.set(VERTEX_COUNT * e.v2 + e.v1);

		}
		cout << segment << endl;
		segments.push_back(segment);
	}
	return segments;
}

void initialize_partition(int* partitioning, vector<int>* partition_members,
		int k) {
	for (int i = 0; i < 4; i++) {
		partitioning[i] = 0;
		partition_members[0].push_back(i);
	}
	for (int i = 4; i < VERTEX_COUNT; i++) {
		partitioning[i] = 1;
		partition_members[1].push_back(i);
	}
}

int* partition_segment(vector<int>* adjacency_list) {
	int* partitioning = new int[VERTEX_COUNT];
	const int k = 2;
	vector<int>* partition_members = new vector<int> [k];

	initialize_partition(partitioning, partition_members, k);

	for (int node = 0; node < VERTEX_COUNT; node++) {
		int best_partition = 0;
		double best_cost = 100000000;
		int i = 0;
		while (partition_members[partitioning[node]][i] != node) {
			i++;
		}
		partition_members[partitioning[node]].erase(partition_members[partitioning[node]].begin()+i);
		for (int partition = 0; partition < k; partition++) {
			partitioning[node] = partition;
			partition_members[partition].push_back(node);
			double cost = node_cost(partitioning, partition_members[partition],
					adjacency_list, node, partition);
			cout << partition << ", " << cost << endl;
			if (cost < best_cost) {
				best_cost = cost;
				best_partition = partition;
			}
			partition_members[partition].erase(
					partition_members[partition].end() - 1);
		}
		partitioning[node] = best_partition;
		partition_members[best_partition].push_back(node);
	}

	return partitioning;
}

int main() {
	map<int, int> nodes;
	map<int, int> renumbered_nodes;
	const string filename = "small_graph.txt";
	vector<int>* adjacency_list = new vector<int> [VERTEX_COUNT];

	vector<int*> partitions;

	priority_queue<Edge, vector<Edge>, CompareEdges> edges = read_graph(
			filename, nodes, renumbered_nodes);

	while (!edges.empty()) {
		Edge e = edges.top();
		while (!edges.empty() && edges.top().timestamp == e.timestamp) {
			e = edges.top();
			edges.pop();
			cout << e.v1 << ", " << e.v2 << ", " << e.timestamp << endl;
			adjacency_list[e.v1].push_back(e.v2);
			adjacency_list[e.v2].push_back(e.v1);

		}
		partitions.push_back(partition_segment(adjacency_list));
	}

	for (int i = 0; i < VERTEX_COUNT; i++) {
		cout << partitions[0][i] << endl;
	}


	return 0;

}
