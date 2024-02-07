#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <algorithm>

using namespace std;

#define df 0.85

vector<vector<size_t>> graph;
vector<int> num_outgoing;
vector<double> pr;
vector<double> new_pr;
vector<double> gather_pr;
vector<double> test;

int start;
int end2;
int num_of_vertex;

template <class Vector, class T>

bool insert_into_vector(Vector& v, const T& t){
    typename Vector::iterator i = lower_bound(v.begin(), v.end(), t);

    if (i == v.end() || t < *i) {
        v.insert(i, t);
        return true;
    } else {
        return false;
    }
}

bool add_arc(size_t from, size_t to){
    vector<size_t> v;
    bool ret = false;
    size_t max_dim = max(from, to);

    if (graph.size() <= max_dim) {
        max_dim = max_dim + 1;
        
        graph.resize(max_dim);
        if (num_outgoing.size() <= max_dim) {
            num_outgoing.resize(max_dim, 0);
        }
    }

    ret = insert_into_vector(graph[to], from);

    if (ret) {
        num_outgoing[from]++;
    }

    return ret;
}

void create_graph_data(string file_path, string delimiter){
    
    istream *infile;

    infile = new ifstream(file_path.c_str());
    size_t line_num = 0;
    string line;
   
   if(infile){
        while(getline(*infile, line)) {
            string from, to;
            size_t pos;

            if(delimiter == " ")
                pos = line.find(" ");
            else if(delimiter == "\t")
                pos = line.find("\t");
            else
                cout << "[ERROR] Input delimiter \" \" or \"t \" to execute program " << endl; 

            from = line.substr(0, pos);
            to = line.substr(pos+1);

            // strtol -> string to int(long)
            add_arc(strtol(from.c_str(), NULL, 10), strtol(to.c_str(), NULL, 10));
            line_num++;
            
            if(line_num % 500000 == 0)
                cerr << "Create " << line_num << " lines" << endl;
      }
   } 
    else {
        cerr << "[ERROR] Unable to open file" <<endl;
        exit(1);
   }

    //add_arc(239822,0);
    num_of_vertex = graph.size();
    cerr << "[ Create " << line_num << " lines, " << num_of_vertex << " vertices graph. ]" << endl;
    
    delete infile;
}

int main(int argc, char* argv[]) {
    double sendbuf;
    double recvbuf_sum;
    struct timespec begin, end1 ;
    struct timespec begin1, end3 ;
    
    // Create graph data;
    string file_path = argv[1];
    create_graph_data(file_path, argv[2]);
 
    cerr << "--------------------------------------------" << endl;
    
    new_pr.resize(num_of_vertex, 0);
    new_pr[0] = 1;
    gather_pr.resize(num_of_vertex, 1/num_of_vertex);

    double dangling_pr = 0;
    double diff = 1;
    double tmp = 0;
    
    double inv_num_of_vertex = 1.0 / num_of_vertex;
    double df_inv = 1.0 - df;

    double* recv_buffer_ptr = gather_pr.data();    
    double* send_buffer_ptr = new_pr.data();

    // vector<vector<size_t>> graph;
    const vector<vector<size_t>>& graph1 = graph;
    const vector<int>& num_outgoing1 = num_outgoing;
    
    cout << "progressing... " << endl;

    clock_gettime(CLOCK_MONOTONIC, &begin1);
    for(int step = 0; step < 10000000; step++){
        tmp = 0;
        
        dangling_pr = 0;
    
        if(step!=0){
            diff = 0;
            for (size_t i = 0; i < num_of_vertex; i++) {
                // outgoing edge가 없으면 
                if (num_outgoing[i] == 0)
                    dangling_pr += gather_pr[i];   
            }
        }
               
        clock_gettime(CLOCK_MONOTONIC, &begin);
        
        for(size_t i = 0; i < num_of_vertex; i++){
            tmp = 0.0;
            
            const size_t graph_size = graph1[i].size();
            const size_t* graph_ptr = graph1[i].data();
            
            for(size_t j=0; j<graph_size; j++){
                const size_t from_page = graph_ptr[j];
                const double inv_num_outgoing = 1.0 / num_outgoing1[from_page];

                tmp += recv_buffer_ptr[from_page]*inv_num_outgoing;
            }
            send_buffer_ptr[i] = (tmp+ dangling_pr*inv_num_of_vertex)*df + df_inv*inv_num_of_vertex;
            diff += fabs(new_pr[i]-gather_pr[i]);
           
        }
        
        clock_gettime(CLOCK_MONOTONIC, &end1);
        long double time1 = (end1.tv_sec - begin.tv_sec) + (end1.tv_nsec - begin.tv_nsec) / 1000000000.0;
        printf("calc 수행시간: %Lfs.\n", time1);

       
        cout << "---------" << step+1 <<"step---------" << endl;
        cout << diff << endl;

        gather_pr = new_pr;
        
        if(diff < 0.00001){
            break;
        }
        
    }
    clock_gettime(CLOCK_MONOTONIC, &end3);
    
    
    
        size_t i;
        double sum = 0;
        cout.precision(numeric_limits<double>::digits10);
         for(i=num_of_vertex - 34;i<num_of_vertex;i++){
            cout << "pr[" <<i<<"]: " << gather_pr[i] <<endl;
        }
        for(i=0;i<num_of_vertex;i++){
            sum += gather_pr[i];
        }
        cout << "s = " <<round(sum) << endl;
        printf("Done.\n");
    
        int important = 0;
        string result = "";
        double important_pr = gather_pr[0];
        double tmp1 = important_pr;
        for (int i=0;i< num_of_vertex;i++){
            important_pr = max(important_pr, gather_pr[i]);
            if(tmp1 != important_pr){
                important = i;
                tmp1 = important_pr;
            }
        }

        cout << "important page is " << important << " and value is " << tmp1 << endl;

        long double time = (end3.tv_sec - begin1.tv_sec) + (end3.tv_nsec - begin1.tv_nsec) / 1000000000.0;
        printf("수행시간: %Lfs.\n", time);


}