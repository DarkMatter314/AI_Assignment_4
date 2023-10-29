#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <unordered_map>
#include <set>
#include<exception>
#include<random>
#include<chrono>
#include <cstdlib> 
#include <cmath>

// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;

// Our graph consists of a list of nodes where each node is represented as follows:
class Graph_Node
{
private:
	string Node_Name;  // Variable name 
	vector<int> Children; // Children of a particular node - these are index of nodes in graph.
	vector<string> Parents; // Parents of a particular node- note these are names of parents
	int nvalues;  // Number of categories a variable represented by this node can take
	vector<string> values; // Categories of possible values
    unordered_map<string,int> value_map; // map from value name to index
	vector<float> CPT; // conditional probability table as a 1-d array . Look for BIF format to understand its meaning

public:
	// Constructor- a node is initialised with its name and its categories
    Graph_Node(string name, int n, vector<string> vals)
	{
		Node_Name=name;
		nvalues=n;
		values=vals;
		for(int i=0;i<vals.size();i++){
			value_map[vals[i]] = i;
		}
	}

	string get_name()
	{
		return Node_Name;
	}

	vector<int> get_children()
	{
		return Children;
	}

	vector<string> get_Parents()
	{
		return Parents;
	}

	vector<float> get_CPT()
	{
		return CPT;
	}

	int get_nvalues()
	{
		return nvalues;
	}

	vector<string> get_values()
	{
		return values;
	}

    int get_value_index(string val_name)
    {
        return value_map[val_name];
    }

	void set_CPT(vector<float> new_CPT)
	{
		CPT.clear();
		CPT=new_CPT;
	}

    void set_Parents(vector<string> Parent_Nodes)
    {
        Parents.clear();
        Parents=Parent_Nodes;
    }

    // add another node in a graph as a child of this node
    int add_child(int new_child_index )
    {
        for(int i=0;i<Children.size();i++)
        {
            if(Children[i]==new_child_index)
                return 0;
        }
        Children.push_back(new_child_index);
        return 1;
    }
};


 // The whole network represted as a list of nodes
class Network{

	vector<Graph_Node> Pres_Graph;
    unordered_map<string,int> node_map;

public:

	int addNode(Graph_Node node)
	{
		Pres_Graph.push_back(node);
        node_map[node.get_name()] = Pres_Graph.size()-1;
		return 0;
	}  

	int netSize()
	{
		return Pres_Graph.size();
	}

    // get the index of node with a given name
    int get_index(string val_name)
    {
        return node_map[val_name];
    }

	// get the node at nth index
    Graph_Node* get_nth_node(int n)
    {
        if (n<Pres_Graph.size()) return &Pres_Graph[n];
        else cerr<<"Error! index greater than size of network\n"; return &Pres_Graph[0];
    }

    //get the iterator of a node with a given name
    int search_node(string val_name)
    {
        for(int i=0;i<Pres_Graph.size();i++){
            if(Pres_Graph[i].get_name().compare(val_name) == 0)
                return i;
        }
    	cout<<"node not found\n";
        return -1;
    }
};

Network read_network(string file_name)
{
	Network Alarm;
	string line;
	int find=0;
  	ifstream myfile(file_name); 
  	string temp;
  	string name;
    vector<string> values;                                                        
  	
    if (myfile.is_open())
    {
    	while (! myfile.eof() )
    	{
    		stringstream ss;
      		getline (myfile,line);
      		ss.str(line);
     		ss>>temp;
     		if(temp.compare("variable")==0)
     		{
     				ss>>name;
     				getline (myfile,line);      
     				stringstream ss2;
     				ss2.str(line);
     				for(int i=0;i<4;i++)
     				{
     					ss2>>temp;	
     				}
     				values.clear();
     				while(temp.compare("};")!=0)
     				{
     					values.push_back(temp);
     					ss2>>temp;
    				}
     				Graph_Node new_node(name,values.size(),values);
     				int pos=Alarm.addNode(new_node);
     		}
     		else if(temp.compare("probability")==0)
     		{ 
				ss>>temp;
				ss>>temp;
				Graph_Node* listIt=Alarm.get_nth_node(Alarm.search_node(temp));
				int index=Alarm.get_index(temp);
				ss>>temp;
				values.clear();
				while(temp.compare(")")!=0)
				{
					Graph_Node* listIt1=Alarm.get_nth_node(Alarm.search_node(temp));
					listIt1->add_child(index);
					values.push_back(temp);
					ss>>temp;
				}
				listIt->set_Parents(values);
				getline (myfile,line);
				stringstream ss2;
				ss2.str(line);
				ss2>>temp;
				ss2>>temp;
				vector<float> curr_CPT;
				string::size_type sz;
				while(temp.compare(";")!=0)
				{
					curr_CPT.push_back(atof(temp.c_str()));	
					ss2>>temp;
				}
				listIt->set_CPT(curr_CPT);
     		}
    	}
    	myfile.close();
  	}
  	return Alarm;
}

double write_output_file(vector<vector<float>> CPT, string file_out_name, string file_in_name)
{
	ifstream infile(file_in_name);
	ofstream outfile(file_out_name);
	string line;
	double total_error = 0; 
	int j = 0;
	cout << "writing to output file now" << endl;
	if(infile.is_open() && outfile.is_open())
	{
		int node_index = 0;
		while(!infile.eof())
		{
			getline(infile, line);
			stringstream ss;
			ss.str(line);
			string temp;
			ss >> temp; //reading the first word.
			if(temp != "probability")
			{
				outfile << line << '\n';
				continue;
			}
			//else we are inside probability.
			outfile << line << '\n'; //writing the probability line.
			getline(infile, line); 
			ss.str(line);
			ss >> temp; //reading the first word.
			while(temp != "table")
			{
				ss >> temp;
			}
			//now we are inside the table.
			outfile << temp << " ";
			for(int i=0; i<CPT[node_index].size(); i++)
			{
				outfile << CPT[node_index][i] << " ";
				ss >> temp; //at the same time we keep getting -1 into the stream.
				total_error += abs(CPT[node_index][i] - atof(temp.c_str()));
			}
			ss >> temp;
			outfile << temp << '\n';
			node_index++; 
		}
	}
	outfile.close();
	infile.close();
	return total_error;
}

class createCPT{

    public:
    string dataFileName;
    Network Alarm;
    int netsize;
    vector<vector<float>> CPT;
	vector<vector<float>> CPT_new; 
	vector<int> total_cnt; 
	vector<int> missing_positions; 
	vector<vector<int>> all_data;
	//vector<float> sample_weights_for_values; //stores the probability of each value of the missing variable.
	unordered_map<string, int> node_name_map; //maps the name of the node to its index in the network.
	vector<bool> markov_blanket_done; //stores the markov blanket of the missing variable.
	vector<set<int>> markov_blanket_set; //stores the markov blanket of the missing variable.
	
	createCPT(Network Alarmcopy, string datafile)
	{
        dataFileName = datafile;
        Alarm = Alarmcopy;
        netsize = Alarm.netSize();
		markov_blanket_done = vector<bool>(netsize, false);
		markov_blanket_set = vector<set<int>>(netsize);
    }

    int CPTindex(vector<int> &cur_data, int index) //TO OPTIMISE. Memoise.
	{
        Graph_Node* currNode = Alarm.get_nth_node(index);
        vector<string> parents = currNode->get_Parents();
        vector<int> parentIndex(parents.size());
        for(int i=0; i<parents.size(); i++){
            parentIndex[i] = Alarm.get_index(parents[i]);
        }
        int cptindex = 0;
        int multiplier = 1;
        for(int i=parentIndex.size()-1; i>=0; i--){
		// cout<<' '<<parents.size()<<" End "<<cptindex<<' ';
			Graph_Node* parentNode = Alarm.get_nth_node(parentIndex[i]);
            if(cur_data[parentIndex[i]] != -1) 
				cptindex += multiplier*(cur_data[parentIndex[i]]);
			
            multiplier *= (Alarm.get_nth_node(parentIndex[i]))->get_nvalues();
        }
		if(cur_data[index] != -1) cptindex += multiplier*(cur_data[index]);
		return cptindex;
        //return ((allVals[index].compare("?") != 0) ? currNode->get_value_index(allVals[index]) : 0) + cptindex*(currNode->get_nvalues());
    }

	float probGivenParents(int index, vector<int> &data)
	{
		int cptindex = CPTindex(data, index);
		Graph_Node* currNode = Alarm.get_nth_node(index);
		int probSamples = CPT[index][cptindex];
		int totalSamples = 0;
		cptindex -= ((CPT[index].size() * data[index]) / (currNode->get_nvalues()));
		for(int i=0; i<currNode->get_nvalues(); i++){
			totalSamples += CPT[index][cptindex + ((CPT[index].size() * i) / (currNode->get_nvalues()))];
		}
		float probability = ((float)probSamples)/((float)totalSamples);
		return probability;
	}

	float probGivenMarkovBlanket(vector<int> &data, int index){ //TO OPTIMISE: use log. - Calculates relative probability
		float probability = probGivenParents(index, data);
		Graph_Node* currNode = Alarm.get_nth_node(index);
		vector<int> children = currNode->get_children();
		for(int i=0; i<children.size(); i++){
			probability *= probGivenParents(children[i], data);
		}
		return probability;
	}

	vector<int> string_to_int_data(vector<string> &allVals) // No data can be missing for this function.
	{
		vector<int> data(netsize);
		for(int i=0; i<netsize; i++){
			data[i] = Alarm.get_nth_node(i)->get_value_index(allVals[i]);
		}
		return data;
	}
	
	vector<float> imputeMissing( int datapoint_index)
	{ 
		int missingIndex = missing_positions[datapoint_index]; //-1 if there is no missing value.
		if(missingIndex == -1){
			cout << "No missing value in this datapoint\n";
			return vector<float>();
		}
		vector<int> data = all_data[datapoint_index];
		if(missingIndex == -1)
		{
			cout << "No missing values in this datapoint\n";
			return vector<float>(1,1);
		}
		float maxProb = 0;
		int maxIndex = 0;
		Graph_Node* missingNode = Alarm.get_nth_node(missingIndex);
		int nvalues = missingNode->get_nvalues();
		vector<float> sampleWeight(nvalues);
		float totalProb = 0;
		for(int i=0; i<nvalues; i++)
		{
			data[missingIndex] = i;
			float currProb = probGivenMarkovBlanket(data, missingIndex);
			sampleWeight[i] = currProb;
			totalProb += currProb;
		}
		for(int i=0; i<nvalues; i++){
			sampleWeight[i] /= totalProb;
		}
		return sampleWeight;
	}

	void store_data()
	{
		ifstream myfile(dataFileName);
		string line;
		int j = 0;
		if(myfile.is_open())
		{
			while(!myfile.eof())
			{
				getline(myfile, line);
				stringstream ss;
				ss.str(line);
				vector<string> vals(netsize);
				all_data.push_back(vector<int>(netsize, -1));
				bool missing = false;
				for(int i=0; i<netsize; i++){
					ss>>vals[i];
					if(vals[i] == "\"?\"")
					{
						missing = true;
						missing_positions.push_back(i); //stores the position of the missing values in the jth line.
					}
					else
					{
						all_data[all_data.size() - 1][i] = (Alarm.get_nth_node(i))->get_value_index(vals[i]); 
					}
				}
				if(!missing)
				{
					missing_positions.push_back(-1); //-1 indicates this position has no missing values.
				}
			}
		}
		cout << "Data stored" << endl; 
	}

	set<int> get_markov_blanket_set(int index)
	{
		if(markov_blanket_done[index]) return markov_blanket_set[index]; //if the markov blanket has already been calculated, return it.
		else
		markov_blanket_done[index] = true; //else,we calculate it and cache it.
		vector<int> markov_blanket;
		Graph_Node* currNode = Alarm.get_nth_node(index);
		vector<string> parent_names = currNode->get_Parents();
		vector<int> parents;
		for(int i=0; i<parent_names.size(); i++){
			parents.push_back(Alarm.get_index(parent_names[i]));
		}
		for(int i=0; i<parents.size(); i++){
			markov_blanket.push_back(parents[i]);
		}
		vector<int> children = currNode->get_children();
		for(int i=0; i<children.size(); i++){
			markov_blanket.push_back(children[i]);
			vector<string> child_parents = Alarm.get_nth_node(children[i])->get_Parents();
			for(int j=0; j<child_parents.size(); j++)
			{
				if(child_parents[j] != currNode->get_name()) //only the other parents of the children.
					markov_blanket.push_back(Alarm.get_index(child_parents[j]));
			}
		}
		markov_blanket_set[index] = set<int>(markov_blanket.begin(), markov_blanket.end()); //we also set this up.
		return markov_blanket_set[index];
	}

	void calculate_probabilities()
	{
		//first we reset our current CPT_new value. 
		for(int i=0; i<netsize; i++)
		{
			for(int j=0; j<CPT_new[i].size(); j++)
			{
				CPT_new[i][j] = 0.01; //initializing as 1 for laplace smoothing.
			}
		}
		
		vector<float> sample_weights_for_values; //stores the probability of each value of the missing variable. Updated for each datapoint.
		for(int datapoint = 0;datapoint < all_data.size(); datapoint++)
		{
			int missing_pos_vals = 0;
			vector<int> cur_data = all_data[datapoint]; //stores the current datapoint.
			// vector<vector<int>> missing_filled_data;
			if(missing_positions[datapoint] >= 0)
			{
				//then there is a missing data in this data entry.
				sample_weights_for_values = imputeMissing(datapoint);
				missing_pos_vals = sample_weights_for_values.size(); 
				for(int i = 0; i < missing_pos_vals; i++)
				{
					cur_data[missing_positions[datapoint]] = i; //updated this value.
					// missing_filled_data.push_back(cur_data); //storing this since we will use it frequently
					int cpt_index = CPTindex(cur_data, missing_positions[datapoint]);
					if(cpt_index >= CPT[missing_positions[datapoint]].size())
					{
						cerr<<"Error! CPTindexVal greater than size of CPT in CPT init()\n";
						// cout<<CPTindexVal<<' '<<CPT[i].size()<<'\n';
						throw exception();
					}
					CPT_new[missing_positions[datapoint]][cpt_index] += sample_weights_for_values[i]; //increasing the count of the CPT by the appropriate value.
				}
				cur_data[missing_positions[datapoint]] = all_data[datapoint][missing_positions[datapoint]]; //resetting the value to its original one for the next steps.
			}

			for(int i=0; i<netsize; i++)
			{
				if(i == missing_positions[datapoint])
				{ 	//Then we do nothing since we have already calculated the probabilities for this in the above block.
					continue;
				}
				else if(get_markov_blanket_set(i).count(missing_positions[datapoint]) == 0) //then the missing value is not in the markov blanket of the current node
				{
					int CPTindexVal = CPTindex(all_data[datapoint], i);
					if(CPTindexVal >= CPT[i].size())
					{
						cerr<<"err not in markov blanket greater than size of CPT\n";
						throw exception();
					}
					CPT_new[i][CPTindexVal] += 1; //Increases the count of the CPT 
					// cout << "missing not in markov blanket\n";
					continue;
				}
				//else, the markov blanket of the current node contains the missing value.
				for(int j = 0; j < missing_pos_vals; j++)
				{
					cur_data[missing_positions[datapoint]] = j; //updated this value.
					int CPTindexVal = CPTindex(cur_data, i);
					if(CPTindexVal >= CPT[i].size())
					{
						cerr<<"CPTindexVal greater than size of CPT\n";
						cout <<"CPTindexVal: " << CPTindexVal << " CPT[i].size(): " << CPT[i].size() << ", and node: " << i << "\n"; 
						cout << "Missing position: " << missing_positions[datapoint] << "\n";
						throw exception(); 
					}
					CPT_new[i][CPTindexVal] += sample_weights_for_values[j];
				}
			}
			// cout << "datapoint " << datapoint << " done";
		}
		CPT = CPT_new; //At the end of the iteration, we update the CPT. 
	}

    void CPTinit()
	{
        CPT.resize(netsize);
		CPT_new.resize(netsize);
		total_cnt.resize(netsize);
        for(int i=0; i<netsize; i++)
		{
            Graph_Node* currNode = Alarm.get_nth_node(i);
            CPT[i].resize(currNode->get_CPT().size()); 
			CPT_new[i].resize(currNode->get_CPT().size());  
			for(int j = 0; j < CPT[i].size(); j++)
			{
				CPT[i][j] = 0.01; //Laplace smoothing. initialising all as 1.
				CPT_new[i][j] = 0.01;
			}   
        }
		this->store_data(); //stores all the data in all_data vector, and the missing values as well.
		int datapoint = 0;
		for(datapoint = 0; datapoint < all_data.size(); datapoint++)
		{
			vector<float> sample_weights_for_values;
			if(missing_positions[datapoint] >= 0)
			{
				//then we must set an initial random value in place of the missing value
				int missing_total_vals = Alarm.get_nth_node(missing_positions[datapoint])->get_nvalues();
				int random_val = rand() % missing_total_vals;
				all_data[datapoint][missing_positions[datapoint]] = random_val; //randomly initializing the missing value.
			}
			for(int i=0; i<netsize; i++)
			{
				int CPTindexVal = CPTindex(all_data[datapoint], i);
				if(CPTindexVal >= CPT[i].size())
				{
					cerr<<"Error! CPTindexVal greater than size of CPT in CPT init()\n";
				}
				CPT[i][CPTindexVal]++; //Increases the count of the CPT 
			}
		}
    }

	void converge_probabilities()
	{
		int iterations = 0;
		while(iterations < 100)
		{
			cout << "Iteration: " << iterations << "\n";
			calculate_probabilities();
			iterations++;
		}
	}

	void convertCountsToCPT()
	{
		vector<vector<float>> finalCPT;
		finalCPT.resize(netsize);
		for(int i=0; i<CPT.size(); i++)
		{
			finalCPT[i].resize(CPT[i].size());
			float total = 0;
			Graph_Node* currNode = Alarm.get_nth_node(i);
			int total_values = currNode->get_nvalues();
			int values_for_one = CPT[i].size() / total_values;
			for(int j=0; j<values_for_one; j++)
			{
				for(int k=0; k < total_values; k++)
				{
					total += CPT[i][k * values_for_one + j];
				}
				for(int k=0; k < total_values; k++)
				{
					finalCPT[i][k * values_for_one + j] = CPT[i][k * values_for_one + j] / total;
				}
			}
		}
		CPT = finalCPT;
	}
};

int main(int argc, char* argv[])
{
	if(argc != 3){
		cout << "Usage: ./em_v1 <bayesnet> <datafile>\n";
		return 0;
	}
	Network Alarm;
	Alarm=read_network(string(argv[1]));
	createCPT CPT(Alarm, string(argv[2]));
	CPT.CPTinit();
	cout<<"Initialised CPT\n";
	CPT.converge_probabilities(); //Converges the probabilities. Currently uses fixed number of 100 iterations.
	CPT.convertCountsToCPT();
	double total_error = write_output_file(CPT.CPT, "solved_counts.bif", string(argv[1]));
	cout << "Total Error: " << total_error << endl;
}
