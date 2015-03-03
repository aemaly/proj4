# proj4
#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cfloat>
#include <unordered_map>
#include <list>
#include <deque>
#include <queue>
#include "getopt.h"

using namespace std;

ostringstream os;

struct Fclty{

    int x_coord; 
    int y_coord;
    int num_fac;
    bool kv;
    double dv;
    int pv;

};

struct Path{

    long long num_fac_start;
    long long num_fac_end;
    double dist;


};

struct TSP{

    double totalweight;
    deque<Fclty> tsplist;

};

class Path_comparator{

public:
    Path_comparator();
    bool operator()(const Path &a, const Path &b) const
    {
        if(a.dist > b.dist) return true;
        else if(a.dist == b.dist && a.num_fac_start > b.num_fac_start) return true;
        else if(a.dist == b.dist && a.num_fac_start == b.num_fac_start
                                      && a.num_fac_end > b.num_fac_end) return true;

        else return false;
    }
};

Path_comparator::Path_comparator() {}

class A{

public:
        unordered_map<int, list<int>> adjacencylist;
        priority_queue<Path, std::vector<Path>,Path_comparator> Pathpicker;
        vector<Fclty> all_Facilities;
        vector<int> unionfinder;
        vector<Path> MST;
        
A() {}

void readA()
{
    string Facilities; int num_facilities = 0;;
    cin>>Facilities>>num_facilities;
    int x = 0; int y = 0;

    Fclty Facility;
    int i;
    for(i = 0; i < num_facilities; ++i)
    {
        cin>>x>>y;
        Facility.x_coord = x; Facility.y_coord = y;
        Facility.num_fac = i; Facility.kv = false;
        Facility.dv = DBL_MAX; Facility.pv = -1;

        all_Facilities.push_back(Facility);
        unionfinder.push_back(i);

    }

    string Paths; int num_paths;
    cin>>Paths>>num_paths;
    int lower_facility_number = 0; int higher_facility_number = 0;

    for(i = 0; i < num_paths; ++ i)
    {
        cin>>lower_facility_number>>higher_facility_number;
        adjacencylist[lower_facility_number].push_back(higher_facility_number);
        //adjacencylist[higher_facility_number].push_back(lower_facility_number);
    }
}

int Find(int x)
{
    if(unionfinder[x] != unionfinder[unionfinder[x]])
        unionfinder[x] = Find(unionfinder[x]);
    return unionfinder[x];
}

bool Union(int start, int end)
{
    int one = Find(start), two = Find(end);
    if(one == two) return false;

    unionfinder[one] = two;
    return true;
}

void KruskMST()
{
    //presort all edges in priority queue
    //loop thru in sorted order
    //if creates cycle, discard
    //initial 2 edges may be disjoint
    //we are growing a forest FROM disjoint trees into one tree
    //ADD EDGE TO TREE

    if(all_Facilities.empty())
    {
        return;
    }
    else
        if(all_Facilities.size() == 1)
        {
            os<<"0"<<endl;
            return;
        }

    //unsigned int i = 0;
    double totaldist = 0;
    unsigned int s = 0;
    for(s = 0; s < adjacencylist.size(); ++s) //unordered, may switch to vector..
    {
        for(auto e = adjacencylist[s].begin(); e != adjacencylist[s].end(); ++e)
        {
            Path newpath;
            newpath.num_fac_start = s; 
            newpath.num_fac_end = *e;
            newpath.dist = sqrt(pow((all_Facilities[s].x_coord -all_Facilities[*e].x_coord),2) +
                                    pow((all_Facilities[s].y_coord - all_Facilities[*e].y_coord),2));
            Pathpicker.push(newpath);
            //cout<<"dist "<<newpath.dist<<endl;
        }
    }


    //use union find algorithm to check for cycles
    unsigned int i = 0;
    while(!(Pathpicker.empty())) //an MST has |v|-1 edges, and allfacilities contains facilities with vertices
    {
        if(Union(Pathpicker.top().num_fac_start, Pathpicker.top().num_fac_end))
        {
            totaldist += Pathpicker.top().dist;
            MST.push_back(Pathpicker.top()); 
            
        }
            Pathpicker.pop();
    }


    os<<totaldist<<endl;
    for(i = 0; i < MST.size(); ++i)
    {
        os<<MST[i].num_fac_start<<" "<<MST[i].num_fac_end<<endl;
    }
}

~A(){}

};

class B{

public:
    vector<Fclty> all_Facilities;

B() {}
    
void readB()
{
    string Facilities; int num_facilities = 0;;
    cin>>Facilities>>num_facilities;
    int x = 0; int y = 0;

    Fclty Facility;
    int i;
    for(i = 0; i < num_facilities; ++i)
    {
        cin>>x>>y;
        Facility.x_coord = x; Facility.y_coord = y;
        Facility.num_fac = i; Facility.kv = false;
        Facility.dv = DBL_MAX; Facility.pv = -1;

        all_Facilities.push_back(Facility);
    }
}

void PrimMST()
{
    if(all_Facilities.empty())
    {
        return;
    }

    else 
        if(all_Facilities.size() == 1)
    {
        os<<"0"<<endl;
        return;
    } 

    else
    {
        all_Facilities[0].dv = 0;
        all_Facilities[0].kv = false;
    }

    Fclty curr_min_fac = all_Facilities[0];
    double new_dist = 0;
    
    unsigned int counter = 0;
    unsigned int i = 0;
    
    double total_dist = 0;

    while(true)
    {
        if(counter == (all_Facilities.size())) //all facilities are now in the in set
            break;

        double curr_min_dist = DBL_MAX;
        unsigned int min_index = 0;
        for(i = counter; i < all_Facilities.size(); ++i)
        {
            if(all_Facilities[i].kv == false && 
                    all_Facilities[i].dv < curr_min_dist)
            {
                min_index = i;
                curr_min_dist = all_Facilities[i].dv;
            }
        }

        total_dist += sqrt(curr_min_dist);
        all_Facilities[min_index].kv = true;

        for(i = counter; i < all_Facilities.size(); ++i)
        {
            if(all_Facilities[i].kv == false)
            {
                new_dist = /*sqrt*/(pow((all_Facilities[min_index].x_coord -all_Facilities[i].x_coord),2) +
                                    pow((all_Facilities[min_index].y_coord - all_Facilities[i].y_coord),2));
                //os<<"i is "<<i<<endl;
                
                if(new_dist < all_Facilities[i].dv)
                {
                    all_Facilities[i].dv = new_dist;
                    all_Facilities[i].pv = all_Facilities[min_index].num_fac;
                }
            }
        }

        //check swapping of all data including parent
        curr_min_fac = all_Facilities[min_index];
        all_Facilities[min_index] = all_Facilities[counter];
        all_Facilities[counter] = curr_min_fac;//make sure
        ++counter;
    }

    os<<total_dist<<endl;
    for(i = 0; i< all_Facilities.size(); ++i)
    {
        if(all_Facilities[i].num_fac < all_Facilities[i].pv)
        {
            os<<all_Facilities[i].num_fac<<" "<<all_Facilities[i].pv<<endl;
        }
        else
        {
            if(all_Facilities[i].pv >= 0)
            os<<all_Facilities[i].pv<<" "<<all_Facilities[i].num_fac<<endl;
        }
    }
    
}

~B(){}

};

class C
{

public:
    //priority_queue<Path, std::vector<Path>,Path_comparator> Pathpicker;
    vector<Fclty> all_Facilities;
    deque<Fclty> vis;
    deque<Fclty> unvis;
    TSP tsp;
    deque<double> weights;
    //unordered_map<int, priority_queue<Path, std::vector<Path>, Path_comparator>> all_Paths;

C() {}

void readC()
{
    string Facilities; long long num_facilities = 0;;
    cin>>Facilities>>num_facilities;
    int x = 0; int y = 0;

    Fclty Facility;
    int i;
    for(i = 0; i < num_facilities; ++i)
    {
        cin>>x>>y;
        Facility.x_coord = x; Facility.y_coord = y;
        Facility.num_fac = i; Facility.kv = false;
        Facility.dv = DBL_MAX; Facility.pv = -1;

        all_Facilities.push_back(Facility);
    }

    tsp.totalweight = DBL_MAX;

    vis.push_back(all_Facilities.front());

    unsigned int j = 1;

    for(j = 1; j < all_Facilities.size(); ++j)
    {
        unvis.push_back(all_Facilities[j]); 
    }
        weights.push_back(0); 
}

void nearestneighbor()
{
    all_Facilities[0].kv = true;
    unsigned int i = 0;
    unsigned int j = 0;
    double weight = 0;
    double nextweight = 0;
    unsigned int curr = 0;
    unsigned int next = 0;
    tsp.tsplist.push_back(all_Facilities[0]);

    for(i = 0; i < all_Facilities.size() -1; ++i)
    {
        double lowest = DBL_MAX;

        for(j = 0; j < all_Facilities.size(); ++j)
        {
            if(all_Facilities[j].kv == false && i != j)
            {
                nextweight = /*sqrt*/(pow((all_Facilities[j].x_coord - all_Facilities[curr].x_coord),2) +
                                    pow((all_Facilities[j].y_coord - all_Facilities[curr].y_coord),2));
            

                if(nextweight < lowest)
                {
                    lowest = nextweight;
                    next = j;
                }
            }

        }

        all_Facilities[next].kv = true;
        curr = next;
        tsp.tsplist.push_back(all_Facilities[curr]);
        weight += sqrt(lowest);

    }
    weight += sqrt(pow((all_Facilities[0].x_coord - all_Facilities[curr].x_coord),2) +
                                    pow((all_Facilities[0].y_coord - all_Facilities[curr].y_coord),2));
    tsp.totalweight = weight;
    
} 

double MSTbound()
{
    double lowerbound = 0;
    if(all_Facilities.empty() || unvis.size() == 0) //useless
    {
        return 0;
    }

    if(unvis.size() == 1)
    {
        return lowerbound;
    }

    else if(unvis.size() == 2)
    {
        lowerbound = sqrt(pow((unvis.front().x_coord - unvis.back().x_coord),2) +
                                    pow((unvis.front().y_coord - unvis.back().y_coord),2)); //length of unvis
    }

    else
    {
        vector<Fclty> new_low(unvis.begin(), unvis.end());
        new_low[0].dv = 0;
        new_low[0].kv = false;
    

    Fclty curr_min_fac = new_low[0];
    double new_dist = 0;
    
    unsigned int counter = 0;
    unsigned int i = 0;
    
    double total_dist = 0;

    while(true)
    {
        if(counter == (new_low.size())) //all facilities are now in the in set
            break;

        double curr_min_dist = DBL_MAX;
        unsigned int min_index = 0;
        for(i = counter; i < new_low.size(); ++i)
        {
            if(new_low[i].kv == false && 
                    new_low[i].dv < curr_min_dist)
            {
                min_index = i;
                curr_min_dist = new_low[i].dv;
            }
        }

        total_dist += sqrt(curr_min_dist);
        new_low[min_index].kv = true;

        for(i = counter; i < new_low.size(); ++i)
        {
            if(new_low[i].kv == false)
            {
                new_dist = /*sqrt*/(pow((new_low[min_index].x_coord -new_low[i].x_coord),2) +
                                    pow((new_low[min_index].y_coord - new_low[i].y_coord),2));
                //os<<"i is "<<i<<endl;
                
                if(new_dist < new_low[i].dv)
                {
                    new_low[i].dv = new_dist;
                    new_low[i].pv = new_low[min_index].num_fac;
                }
            }
        }
        //check all data including parent
        curr_min_fac = new_low[min_index];
        new_low[min_index] = new_low[counter];
        new_low[counter] = curr_min_fac;//make sure
        ++counter;
    }

    lowerbound = total_dist; //use lowbound as lower bound for TSP
    }
//connect the bounds in permute function
    return lowerbound;
}

void permute()
{
    double nextweight = 0;
    double realbound = 0;
    unsigned int i = 0;
    if(unvis.empty()) //complete path!
    {
        nextweight = weights.front() +
            sqrt(pow((vis.front().x_coord -unvis.back().x_coord),2) +
                                    pow((vis.front().y_coord - unvis.back().y_coord),2)); //should make function lol

    //cout<<"hi"<<endl;
        if (nextweight < tsp.totalweight)
        {
            tsp.tsplist = vis;
            tsp.totalweight = nextweight;
        }

        return;
    }

    for(i = 0; i < unvis.size(); ++i)
    {
        nextweight = weights.front() +
            sqrt(pow((vis.back().x_coord -unvis.front().x_coord),2) +
                                    pow((vis.back().y_coord - unvis.front().y_coord),2));
        weights.push_front(nextweight); //update weight

        vis.push_back(unvis.front()); //mark as visited
        unvis.pop_front(); //now visited

        realbound = weights.front() + MSTbound(); //finds unvisited lowerbound
        if(realbound < tsp.totalweight)
            permute();

        unvis.push_back(vis.back());
        vis.pop_back();
        weights.pop_front();
    }
}

void PrimTSP()
{

    if(all_Facilities.empty() || all_Facilities.size() == 1 || all_Facilities.size() == 2)
        return;

        nearestneighbor();
        permute();

    os<<tsp.totalweight<<endl;
    unsigned int i = 0;
    for(i = 0; i < tsp.tsplist.size(); ++i)
        os<< tsp.tsplist[i].num_fac<<" ";
    os<<endl;
}

~C(){}

};

int main(int argc, char** argv)
{
    ios_base::sync_with_stdio();
    os<<std::setprecision(2);
    os<<std::fixed;

struct option long0pts[] = {

    { "clientType", required_argument, NULL, 'c' },
    { "help", no_argument, NULL, 'h' }
                            };

    opterr = false;
    bool client = false;
    char* CLIENT_TYPE = nullptr;

    int opt = 0, index = 0;

    while((opt = getopt_long(argc, argv, "c:h", long0pts, &index)) != -1)
    {
        switch (opt) 
        {
        case 'c': 
            client = true;
            CLIENT_TYPE = optarg;
            break;
        case 'h':
            cerr<<"Minimum spanning tree. That is all."<<endl;
            exit(0);
            break;

        default:
            cerr<< "I didn't recognize one of your flags!"<<endl;
            exit(1);
        }
    }

        if(client == false)
        {
            cerr<<"No client or type specified"<<endl;
            exit(1);
        }

        else if(!(strcmp(CLIENT_TYPE, "A") == 0) && 
                    !(strcmp(CLIENT_TYPE, "B") == 0) && 
                        !(strcmp(CLIENT_TYPE, "C") == 0)) //error check getopt
        {
            cerr<<"You didn't specify a valid client type"<<endl;
            exit(1);
        }


    if(strcmp(CLIENT_TYPE, "A") == 0)
    {
        A a;
        a.readA();
        a.KruskMST();
    }

    else if(strcmp(CLIENT_TYPE,"B") == 0)
    {
        B b;
        b.readB();
        b.PrimMST();
    }

    else
    {
        C c;
        c.readC();
        c.PrimTSP();
    }

cout<<os.str();
return 0;
}
