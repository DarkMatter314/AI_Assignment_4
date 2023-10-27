#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <unordered_map>


// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;

// Our graph consists of a list of nodes where each node is represented as follows:
class Graph_Node{

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
    Graph_Node(string name,int n,vector<string> vals)
	{
		Node_Name=name;
		nvalues=n;
		values=vals;
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
        else cout<<"Error! index greater than size of network\n"; return &Pres_Graph[0];
    }

    //get the iterator of a node with a given name
    int search_node(string val_name)
    {
        for(int i=0;i<Pres_Graph.size();i++){
            if(Pres_Graph[i].get_name().compare(val_name)==0)
                return i;
        }
    	cout<<"node not found\n";
        return -1;
    }
};

Network read_network()
{
	Network Alarm;
	string line;
	int find=0;
  	ifstream myfile("./data/alarm.bif"); 
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

class createCPT{

    public:
    string dataFileName;
    Network Alarm;
    int netsize;
    vector<float> CPT;

	createCPT(Network Alarmcopy, string datafile){
        dataFileName = datafile;
        Alarm = Alarmcopy;
        netsize = Alarm.netSize();
    }

    int CPTindex(vector<string> &allVals, int index){
        Graph_Node* currNode = Alarm.get_nth_node(index);
        vector<string> parents = (Alarm.get_nth_node(index))->get_Parents();
        vector<int> parentIndex(parents.size());
        for(int i=0; i<parents.size(); i++){
            parentIndex[i] = Alarm.get_index(parents[i]);
        }
        int cptindex = 0;
        int multiplier = 1;
        for(int i=parentIndex.size()-1; i>=0; i--){
            if(!allVals[parentIndex[i]].compare("?")) cptindex += multiplier*(currNode->get_value_index(allVals[parentIndex[i]]));
            multiplier *= (Alarm.get_nth_node(parentIndex[i]))->get_nvalues();
        }
        return currNode->get_value_index(allVals[index]) + cptindex*(currNode->get_nvalues());
    }

    void CPTinit(){
        vector<vector<float>> CPT(netsize);
        for(int i=0; i<netsize; i++){
            Graph_Node* currNode = Alarm.get_nth_node(i);
            CPT[i].resize(currNode->get_CPT().size());      
        }
        ifstream myfile(dataFileName);
        string line;
        if(myfile.is_open()){
            while(!myfile.eof()){
                getline(myfile, line);
                stringstream ss;
                ss.str(line);
                vector<string> allVals(netsize);
                for(int i=0; i<netsize; i++){
                    ss>>allVals[i];
                    allVals[i] = allVals[i].substr(1, allVals[i].size()-2);
                }
                for(int i=0; i<netsize; i++){
                    CPT[i][CPTindex(allVals, i)]++;
                }
            }
        }
    }

};

int main()
{
	Network Alarm;
	Alarm=read_network();
	// Example: to do something
	cout<<"Perfect! Hurrah! \n";
}




