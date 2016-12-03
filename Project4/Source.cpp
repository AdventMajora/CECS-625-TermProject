#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <time.h>

using namespace std;

struct city {
	int city_name;
	float x;
	float y;

	city() {

	}

	city(int name, float x_coord, float y_coord) {
		city_name = name;
		x = x_coord;
		y = y_coord;
	}
};

struct path {
	vector<city> route;
	float distance;
	float prob;

	path() {

	}

	path(vector<city> new_route, float route_distance) {
		route = new_route;
		distance = route_distance;
	}

	bool operator < (const path& compared_path) const
	{
		return (distance < compared_path.distance);
	}
};

vector<string> read_file(string);
vector<string> split(string&, char);
void gen_dist_table();
float get_distance(city, city);
float route_distance(vector<city>);
path gen_path();
void gen_populations();
void sim_generation();
vector<city> crossover(int, int, int);
vector<city> mutation(vector<city>);

int num_pops = 1;
int pop_size = 100;
int dimension;
string edge_weight;
vector<city> cities;
vector<vector<float>> dist_table;
vector<vector<path>> populations;
int current_generation = 0;
int generation_limit = 500;
string test_file = "Random22.tsp";
int mutation_prob = 90;
int elitism = .05;

int main() {
	cout << "Hellow orld!" << endl;
	srand(time(NULL));

	//read in scenario
	vector<string> settings = read_file(test_file);

	//parse input for simulation settings
	for (unsigned i = 0; i < settings.size(); i++) {
		string line = settings[i];
		if ((int)line.find("DIMENSION") > -1) {
			dimension = stoi(line.substr(line.find("DIMENSION") + 11));
		}
		if ((int)line.find("EDGE") > -1) {
			edge_weight = line.substr(line.find("EDGE") + 18);
		}
		if ((int)line.find("NODE") > -1) {
			vector<string> details;
			for (unsigned j = i + 1; j < i + dimension + 1; j++) {
				details = split(settings[j], ' ');
				cities.push_back(city(stoi(details[0]),stof(details[1]), stof(details[2])));
			}
		}
	}
	cout << "Dimension: " << dimension << endl;
	cout << "Edge Weight: " << edge_weight << endl;
	cout << "Cities: " << cities.size() << endl;
	
	//generate distance table
	gen_dist_table();
	
	//generate population
	gen_populations();
	
	/*cout << "Initial population: " << endl;
	for (unsigned i = 0; i < populations[0].size(); i++) {
		for (unsigned j = 0; j < populations[0][i].route.size(); j++) {
			cout << populations[0][i].route[j].city_name << ", ";
		}
		cout << populations[0][i].distance << endl;
	}*/

	//do the things
	while (current_generation < generation_limit) {
		sim_generation();
	}

	/*cout << "resultant population: " << endl;
	for (unsigned i = 0; i < populations[0].size(); i++) {
		for (unsigned j = 0; j < populations[0][i].route.size(); j++) {
			cout << populations[0][i].route[j].city_name << ", ";
		}
		cout << populations[0][i].distance << endl;
	}*/

	cout << "best result: " << endl;
	for (unsigned j = 0; j < populations[0][0].route.size(); j++) {
		cout << populations[0][0].route[j].city_name << ", ";
	}
	cout << populations[0][0].distance << endl;

	system("pause");
	return 0;
}

//read in a file
vector<string> read_file(string filepath) {
	vector<string> contents;
	string line;
	ifstream myfile(filepath);
	if (myfile) {
		while (getline(myfile, line)) {
			contents.push_back(line.c_str());
		}
		myfile.close();
	} else {
		cout << filepath << "was unable to be opened..." << endl;
	}
	return contents;
}

//split a strigh on a delimiter
vector<string> split(string &s, char delim) {
	vector<string> elems;
	stringstream ss;
	ss.str(s);
	string item;
	while (getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

//generate a lookup table for distances between cities
void gen_dist_table() {
	for (int i = 0; i < dimension; i++) {
		vector<float> col;
		for (int j = 0; j < dimension; j++) {
			if (j == i) {
				col.push_back(-1.0);
			} else {
				col.push_back(get_distance(cities[i], cities[j]));
			}
		}
		dist_table.push_back(col);
	}
}

//calculate the distance between two cities
float get_distance(city a, city b) {
	return sqrt(abs(pow((b.x-a.x),2)-pow((b.y-a.y),2)));
}

//generate the initial population(s)
void gen_populations() {
	//generate a population
	for (int i = 0; i < num_pops; i++) {
		vector<path> new_pop;
		path new_path;
		for (int j = 0; j < pop_size; j++) {
			new_path = gen_path();
			new_pop.push_back(new_path);
		}
		//sort it in ascending order
		sort(new_pop.begin(), new_pop.end());

		populations.push_back(new_pop);
	}
}

//generates a rendom valid path
path gen_path() {
	path new_path;
	vector<city> to_shuffle = cities;
	int m = to_shuffle.size();
	int i;
	city t;
	while (m > 0) {
		m--;
		i = (int)floor((((double)rand()/(RAND_MAX))) *m);
		t = to_shuffle[m];
		to_shuffle[m] = to_shuffle[i];
		to_shuffle[i] = t;
	}
	to_shuffle.push_back(to_shuffle[0]);
	new_path.route = to_shuffle;
	new_path.distance = route_distance(to_shuffle);
	return new_path;
}

//calculates the total distance along a path
float route_distance(vector<city> calc_path) {
	float total_distance = 0;
	for (unsigned i = 0; i < calc_path.size()-1; i++) {
		float branch = dist_table[calc_path[i].city_name - 1][calc_path[i + 1].city_name - 1];
		total_distance += branch;
	}
	return total_distance;
}

//simulates a generation of growth
void sim_generation() {
	if (current_generation > generation_limit) {
		cout << "Simultion complete!" << endl;
	} else {
		float total_distances = 0;
		float running_distance = 0;
		float total_probs = 0;
		int p1_path_index, p2_path_index;
		int sel;
		float rnd;
		string tmp;
		vector<city> offspring;
		vector<path> new_pop, final_pop;

		current_generation++;
		for (unsigned p = 0; p < populations.size(); p++) {
			new_pop.clear();
			//PARALLELIZE THIS SHIT HERE
			for (unsigned i = 0; i < populations[p].size(); i++) {
				
				//set up roulette
				total_distances = 0;
				running_distance = 0;
				total_probs = 0;
				for (unsigned j = 0; j < populations[p].size(); j++) {
					total_distances += populations[p][j].distance;
				}
				for (unsigned j = 0; j < populations[p].size(); j++) {
					populations[p][j].prob = 1 - (populations[p][j].distance / total_distances);
					total_probs += populations[p][j].prob;
				}

				//select parent 1
				rnd = (int)floor(((double)rand() / (RAND_MAX))*total_probs);
				for (sel = 0; sel < (int)populations[p].size() && rnd > 0; sel++) {
					rnd -= populations[p][sel].prob;
				}
				if (sel < 1) {
					sel=1;
				}
				p1_path_index = sel-1;

				//select parent 2
				rnd = floor(((double)rand() / (RAND_MAX))*total_probs);
				for (sel = 0; sel < (int)populations[p].size() && rnd > 0; sel++) {
					rnd -= populations[p][sel].prob;
				}
				if (sel < 1) {
					sel=1;
				}
				p2_path_index = sel - 1;

				offspring.clear();

				//perform crossover
				offspring = crossover(p1_path_index, p2_path_index, p);

				//roll for chance to mutate
				if (floor(((double)rand() / (RAND_MAX)) * 100 + 1) < mutation_prob) {
					offspring = mutation(offspring);
				}

				new_pop.push_back(path(offspring, route_distance(offspring)));
			}
			sort(new_pop.begin(), new_pop.end());

			//hey ho! it's elitism!
			final_pop.clear();
			int elites = (int)(pop_size*.1);
			if (elites < 1) {
				elites = 1;
			}
			for (int i = 0; i < elites; i++) {
				final_pop.push_back(populations[p][i]);
			}
			for (unsigned i = 0; i < new_pop.size()-elites; i++) {
				final_pop.push_back(new_pop[i]);
			}

			sort(final_pop.begin(), final_pop.end());

			populations[p] = final_pop;
			cout << populations[0][0].distance << endl;
		}
	}
}

//performs genetic crossover
vector<city> crossover(int a, int b, int p) {
	vector<city> offspring;
	int tmp;

	//a[5,2,4,3,6,1,5]
	//b[1,2,3,4,5,6,1]
	//o[]

	int cutoff = (int)floor((((double)rand() / (RAND_MAX))*populations[p][a].route.size() - 2));
	if (cutoff > populations[p][a].route.size() - 2) {
		cutoff = populations[p][a].route.size() - 2;
	}
	if (cutoff < 1) {
		cutoff = 1;
	}
	
	//build upto the cutoff
	for (int j = 0; j < cutoff; j++) {
		offspring.push_back(populations[p][a].route[j]);
	}

	//a[5,2,4,3,6,1,5]
	//b[1,2,3,4,5,6,1]
	//o[5,2,4]

	//make a "copy" of offspring for comparisons later
	vector<int> offstring;
	for (unsigned j = 0; j < offspring.size(); j++) {
		offstring.push_back((offspring[j].city_name));
	}

	//make a "copy" of b's path
	vector<int> p2_copy;
	for (unsigned j = 0; j < populations[p][b].route.size(); j++) {
		p2_copy.push_back((populations[p][b].route[j].city_name));
	}

	//re-order repeated nodes
	for (int j = 0; j < cutoff; j++) {
		if (p2_copy[j] != offspring[j].city_name) {
			int oldLoc = find(p2_copy.begin(), p2_copy.end(), offspring[j].city_name) - p2_copy.begin();
			tmp = p2_copy[j];
			p2_copy[j] = p2_copy[oldLoc];
			p2_copy[oldLoc] = tmp;
		}
	}

	//b[5,2,4,3,1,6,1]
	//o[5,2,4]

	//get the rest of the path from parent 2
	for (unsigned j = cutoff; j < p2_copy.size() - 1; j++) {
		offspring.push_back(cities[(p2_copy[j])-1]);
	}

	//b[5,2,4,3,1,6,1]
	//o[5,2,4,3.1.6]

	//loop the path
	offspring.push_back(offspring[0]);

	return offspring;
}

//performs genetic mutation
vector<city> mutation(vector<city> route) {

	vector<city> mutated = route;
	int mut_start = (int)floor(((double)rand() / (RAND_MAX))*(mutated.size() - 2)) + 1;
	int mut_end = (int)floor(((double)rand() / (RAND_MAX))*(mutated.size() - 2)) + mut_start;

	if (mut_end > mutated.size() - 2) {
		mut_end = mutated.size() - 3;
	}

	vector<city> sub;
	for (int i = mut_end-1; i >= mut_start; i--) {
		sub.push_back(mutated[i]);
	}

	for (int i = mut_start; i < mut_end; i++) {
		mutated[i] = sub[i - mut_start];
	}
	
	return mutated;
}
