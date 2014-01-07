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

int VERTEX_COUNT;
vector<int>* adjacency_list;
int* partitioning;
vector<vector<int> > partition_members;

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
	int* partition;
	int timestamp;
	GraphSegment(int* partition, int timestamp) :
			partition(partition), timestamp(timestamp) {
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

priority_queue<Edge, vector<Edge>, CompareEdges> read_graph(string filename) {
	priority_queue<Edge, vector<Edge>, CompareEdges> edges;
	int max_node_index = 0;
	string line;
	ifstream myfile(filename.c_str());
	if (myfile.is_open()) {
		while (getline(myfile, line)) {
			vector<string> splitRow = split(line, ';');
			int v1 = atoi(splitRow[0].c_str());
			int v2 = atoi(splitRow[1].c_str());
			int timestamp = atoi(splitRow[2].c_str());
			edges.push(Edge(v1, v2, timestamp));
			if (v1 > max_node_index) {
				max_node_index = v1;
			}
			if (v2 > max_node_index) {
				max_node_index = v2;
			}
		}
	}
	myfile.close();
	VERTEX_COUNT = max_node_index + 1;
	return edges;
}

double node_cost(int* partitioning, vector<vector<int> > partition_members,
		int node, int partition) {
	const int k = partition_members.size();

	int possible_edges_from_node = VERTEX_COUNT - 1;
	int partition_size = partition_members[partition].size();
	int possible_edges_from_partition = partition_size
			* (VERTEX_COUNT - partition_size)
			+ partition_size * (partition_size - 1);
	int node_neighbors_count = adjacency_list[node].size();
	int partition_neighbors_count = 0;
	int partition_neighbors_in_same_partition = 0;

	double* node_neighbors_in_partition = new double[k];
	double* partition_neighbors_in_partition = new double[k];
	for (int i = 0; i < k; i++) {
		node_neighbors_in_partition[i] = 0;
		partition_neighbors_in_partition[i] = 0;
	}

	for (vector<int>::iterator it = adjacency_list[node].begin();
			it != adjacency_list[node].end(); ++it) {
		int neighbor = *it;
		node_neighbors_in_partition[partitioning[neighbor]]++;
	}

	for (vector<int>::iterator it = partition_members[partition].begin();
			it != partition_members[partition].end(); it++) {
		for (vector<int>::iterator it2 = adjacency_list[*it].begin();
				it2 != adjacency_list[*it].end(); it2++) {
			int neighbor = *it2;
			if (partitioning[neighbor] == partition) {
				partition_neighbors_in_same_partition++;
			} else {
				partition_neighbors_in_partition[partitioning[neighbor]]++;
				partition_neighbors_count++;
			}
		}
	}
	partition_neighbors_in_partition[partition] +=
			partition_neighbors_in_same_partition / 2;
	partition_neighbors_count += partition_neighbors_count / 2;
//	cout << "node " << node + 1 << ", partition " << partition + 1 << endl;
//	for (int i = 0; i < k; i++) {
//		cout << "partition " << i + 1 << ": " << node_neighbors_in_partition[i]
//				<< ", " << partition_neighbors_in_partition[i] << endl;
//	}

	double entropy = 0;
	for (int i = 0; i < k; i++) {
		if (node_neighbors_in_partition[i] > 0) {
			entropy += (node_neighbors_in_partition[i]
					/ possible_edges_from_node)
					* (log(
							partition_neighbors_in_partition[i]
									/ possible_edges_from_partition) / log(2));
		}
	}
	if (possible_edges_from_partition != partition_neighbors_count) {
		entropy += ((double) (possible_edges_from_node - node_neighbors_count)
				/ possible_edges_from_node)
				* (log(
						(double) (possible_edges_from_partition
								- partition_neighbors_count)
								/ possible_edges_from_partition) / log(2));
	}
	delete node_neighbors_in_partition;
	delete partition_neighbors_in_partition;
	return -entropy;
}

double entropy(int* neighbors_in_partition, int k) {
	int possible_edges_from_node = VERTEX_COUNT - 1;
	double entropy = 0;
	int edges_count = 0;
	for (int i = 0; i < k; i++) {
		if (neighbors_in_partition[i] > 0) {
			entropy += ((double) neighbors_in_partition[i]
					/ (double) possible_edges_from_node)
					* (log(
							(double) neighbors_in_partition[i]
									/ (double) possible_edges_from_node)
							/ log(2));
			edges_count += neighbors_in_partition[i];
		}
	}
	if (possible_edges_from_node != edges_count) {
		entropy += ((double) (possible_edges_from_node - edges_count)
				/ possible_edges_from_node)
				* (log(
						(double) (possible_edges_from_node - edges_count)
								/ possible_edges_from_node) / log(2));
	}
	return -entropy;
}

double average_entropy(int* partitioning,
		vector<vector<int> > partition_members, int partition) {
	double partition_entropy = 0;
	for (vector<int>::iterator it = partition_members[partition].begin();
			it != partition_members[partition].end(); it++) {
		int* neighbors_in_partition = new int[partition_members.size()];
		for (unsigned int i = 0; i < partition_members.size(); i++) {
			neighbors_in_partition[i] = 0;
		}
		for (vector<int>::iterator iter = adjacency_list[*it].begin();
				iter != adjacency_list[*it].end(); iter++) {
			neighbors_in_partition[partitioning[*iter]]++;
		}
		partition_entropy += entropy(neighbors_in_partition,
				partition_members.size());
	}
	return (partition_entropy / partition_members[partition].size());
}

int find_partition_with_largest_entropy() {
	double smallest_average_entropy = 100000000;
	int best_partition = -1;
	for (unsigned int partition = 0; partition < partition_members.size();
			partition++) {
		double partition_average_entropy = average_entropy(partitioning,
				partition_members, partition);
		if (partition_average_entropy < smallest_average_entropy) {
			smallest_average_entropy = partition_average_entropy;
			best_partition = partition;
		}
	}
	return best_partition;
}

/*
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
 */

void initialize_partition() {
	partitioning = new int[VERTEX_COUNT];
	partition_members = vector<vector<int> >(1);
	for (int i = 0; i < VERTEX_COUNT; i++) {
		partitioning[i] = 0;
		partition_members[0].push_back(i);
	}
}

double total_cost(vector<int>* adjacency_list, int* partitioning,
		vector<vector<int> > partition_members, const int k) {
	double cost = 0;
	for (int i = 0; i < k; i++) {
		for (vector<int>::iterator it = partition_members[i].begin();
				it != partition_members[i].end(); it++) {
			int node = *it;
			cost += node_cost(partitioning, partition_members, node,
					partitioning[node]);
		}

	}
	return cost;
}

bool update_partitions() {
	const int k = partition_members.size();
	bool changed = false;
	for (int node = 0; node < VERTEX_COUNT; node++) {
		if (partition_members[partitioning[node]].size() > 1) {
			int best_partition = 0;
			double best_cost = 100000000;
			int i = 0;
			while (partition_members[partitioning[node]][i] != node) {
				i++;
			}
			partition_members[partitioning[node]].erase(
					partition_members[partitioning[node]].begin() + i);
			for (int partition = 0; partition < k; partition++) {
				partitioning[node] = partition;
				partition_members[partition].push_back(node);
				double cost = node_cost(partitioning, partition_members, node,
						partition);
				if (cost < best_cost) {
					best_cost = cost;
					best_partition = partition;
				}
				partition_members[partition].erase(
						partition_members[partition].end() - 1);
			}
			if (best_partition != partitioning[node]) {
				changed = true;
			}
			partitioning[node] = best_partition;
			partition_members[best_partition].push_back(node);
		}
	}

	return changed;
}

int* searchK() {
	// TODO: now only one iteration of split, update, merge. Need to make repeatable, but it stayed in infinite loop:
	// merge merged into [(4, 6, 7, 5, )(1, 2, 3, )] and split split into [(4, 6, 7, )(1, 2, 3, )(5, )]
	bool changed = true;
	int k = partition_members.size();
//	while (changed) {
	changed = false;

	// try to split
	int partition = find_partition_with_largest_entropy();
	double current_average_entropy = average_entropy(partitioning,
			partition_members, partition);
	unsigned int position = 0;
	while (position < partition_members[partition].size()) {
		int node = partition_members[partition][position];
		partitioning[node] = k;
		vector<vector<int> > new_partition_members = partition_members;
		new_partition_members[partition].erase(
				new_partition_members[partition].begin() + position);
		new_partition_members.push_back(vector<int>());

		new_partition_members[k].push_back(node);

		double new_average_entropy = average_entropy(partitioning,
				new_partition_members, partition);

		if (current_average_entropy - new_average_entropy > 0.0001) {
			k++;
			changed = true;
			partition_members = new_partition_members;
			current_average_entropy = new_average_entropy;
		} else {
			partitioning[node] = partition;
			position++;

		}
	}
	cout << "after split: [";
	for (int j = 0; j < k; j++) {
		cout << "(";
		for (vector<int>::iterator it = partition_members[j].begin();
				it != partition_members[j].end(); it++) {
			cout << (*it) << ", ";
		}
		cout << ")";
	}
	cout << "]" << endl;

	// update partitions
	changed = changed | update_partitions();
	cout << "after update: [";
	for (int j = 0; j < k; j++) {
		cout << "(";
		for (vector<int>::iterator it = partition_members[j].begin();
				it != partition_members[j].end(); it++) {
			cout << (*it) << ", ";
		}
		cout << ")";
	}
	cout << "]" << endl;

	// try to merge
	double current_total_cost = total_cost(adjacency_list, partitioning,
			partition_members, k);
	for (int partition1 = 0; partition1 < k; partition1++) {
		for (int partition2 = partition1 + 1; partition2 < k; partition2++) {
			int* new_partitioning = new int[VERTEX_COUNT];
			for (int node = 0; node < VERTEX_COUNT; node++) {
				if (partitioning[node] == partition2) {
					new_partitioning[node] = partition1;
				} else if (partitioning[node] > partition2) {
					new_partitioning[node] = partitioning[node] - 1;
				} else {
					new_partitioning[node] = partitioning[node];
				}
			}
			vector<vector<int> > new_partition_members = partition_members;
			new_partition_members[partition1].insert(
					new_partition_members[partition1].end(),
					new_partition_members[partition2].begin(),
					new_partition_members[partition2].end());

			new_partition_members.erase(
					new_partition_members.begin() + partition2);
			double new_total_cost = total_cost(adjacency_list, new_partitioning,
					new_partition_members, k - 1);

			if (current_total_cost > new_total_cost) {
				partitioning = new_partitioning;
				partition_members = new_partition_members;
				k = k - 1;
				changed = true;
				current_total_cost = new_total_cost;
				partition2--;
			}
		}
	}
	cout << "after merge: [";
	for (int j = 0; j < k; j++) {
		cout << "(";
		for (vector<int>::iterator it = partition_members[j].begin();
				it != partition_members[j].end(); it++) {
			cout << (*it) << ", ";
		}
		cout << ")";
	}
	cout << "]" << endl;

	return partitioning;
}

int* graphScope() {
	return searchK();
}

int main() {
	cout
			<< "*note that two different representations of partitionings are in use:"
			<< endl;
	cout
			<< "[(3, 5, 6, 4, )(0, 1, 2, )] shows that there are cluster 0 with nodes 3, 5, 6 and 4 and cluster 1 with nodes 0, 1 and 2"
			<< endl;
	cout
			<< "[1, 1, 1, 0, 0, 0, 0] shows that nodes 0, 1 and 2 belong to the 1st cluster and 3, 4, 5 and 6 belong to 0th cluster"
			<< endl;
	cout << endl;

	const string filename = "extractedNetwork.o43237";

	// read network and initialize partitioning (one big cluster)
	priority_queue<Edge, vector<Edge>, CompareEdges> edges = read_graph(
			filename);
	adjacency_list = new vector<int> [VERTEX_COUNT];
	initialize_partition();
	vector<GraphSegment> segments;

	// add edges to adjacency list and find partitioning for each different timestep
	while (!edges.empty()) {
		Edge e = edges.top();
		while (!edges.empty() && edges.top().timestamp == e.timestamp) {
			e = edges.top();
			edges.pop();
			cout << e.v1 << ", " << e.v2 << ", " << e.timestamp << endl;
			adjacency_list[e.v1].push_back(e.v2);
			adjacency_list[e.v2].push_back(e.v1);
		}
		int* segment = graphScope();
		cout << "[";
		for (int i = 0; i < VERTEX_COUNT; i++) {
			cout << segment[i];
			if (i < VERTEX_COUNT - 1) {
				cout << ", ";
			}
		}
		cout << "]" << endl;
		segments.push_back(GraphSegment(segment, e.timestamp));
		cout << endl;
	}

	// print partitionings for each timestamp
	for (vector<GraphSegment>::iterator iter = segments.begin();
			iter != segments.end(); iter++) {
		GraphSegment segment = *iter;
		cout << "time " << segment.timestamp << ": [";
		for (int i = 0; i < VERTEX_COUNT; i++) {
			cout << (segment.partition)[i];
			if (i < VERTEX_COUNT - 1) {
				cout << ", ";
			}
		}
		cout << "]" << endl;
	}

	return 0;

}
