/*
ALgorithm FME
	Repeat
		choose an unkown xj to project
		L = set of row indices i, where aij <0
		U = set of row indices i, where aij >0
		if L or U is empty, its an unconstrained variable.. exit
		for i that belongs to L union U do
			divide the whole row Ai and element bi by absolute value of Aij ( where j is the chosen varaible to eliminate)
		end for
		for i belongs to L do
			for k belongs to U do
				add a new inequality row Ai + Row Ak <= bi +bk
			end for
		end for
		ignore all rows in L and U for further consideration
							
	until done
end

*/
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <cstdlib>
#include <cmath>

using namespace std;
double get_min(vector<double> v  ){
    // returns minimum value of a vector
    double min = v[0];
    for (int i=0; i<v.size(); i++){
        if (v[i]<min) min = v[i];
    }
    return min;
}
double get_max(vector<double> v  ){
    // returns maximum value of a vector
    double max = v[0];
    for (int i=0; i<v.size(); i++){
        if (v[i]>max) max = v[i];
    }
    return max;
}

void printlist(vector<double>  &A){
cout<<"Print list"<<endl;
		for(auto j: A){
			cout<<j<<" ";
		}
		cout<<endl;
	cout<<"Print list end"<<endl;
	return;
}

void printmatrix(vector<vector<double> > &A){
cout<<"Print Matrix"<<endl;
	for(size_t i=0; i<A.size(); i++){
		for(auto j: A[i]){
			cout<<j<<" ";
		}
		cout<<endl;
	}
	cout<<"Print Matrix end"<<endl;
	return;
}
void printset(set<int>  &s)
{ cout<<"Print set"<<endl;
    for (auto  &i: s) {
        cout << i << " ";
    }cout<<"Print set end"<<endl;
    return;
}
vector<vector<double>> FME( int rows, int columns, vector<vector<double>> A,vector<double> b){
			
	// j is variable we choose to eliminate, let j = columns -1, we start eliminating from the last variable (targetting inner most loop)
	int j = columns-1,i,k,lb =0,ub =0;
	vector<vector<double>>  sol ;
	set<int> L,U,Z;
	set<int>::iterator itr,itr1;
	vector<vector<double>> newA ;
	vector<double> newb,tempRowA;
	// FInd the lower & upper bounds for the column j
	for(int i=0;i<rows;i++){
		// if Aij < 0 the index of this row is considered as lowerbound
		if(A[i][j] < 0){
			L.insert(i);
		}
		// if Aij > 0 the index of this row is considered as upperbound
		else if (A[i][j] > 0){
			U.insert(i);
		} else {
			// when the value is zero, we will not consider the row for normalization but add it later
			// some implementations recommend to discard values like this, but we wont be doing it for now.
			Z.insert(i);
		}
	}
//	cout<<"L \n";
//	printset(L);
//	cout<<" U \n";
//	printset(U);
//	cout<<"Z \n";
//	printset(Z);
	
	//if L or U is empty, its an unconstrained variable.. exit
	if(L.empty() || U.empty()){
		cout<<"unconstrained variable x"<<j+1<<endl<<"No closed iteration space"<<endl;
		return sol;
		
	}
	//Proceed to normalization
	if( !L.empty() ){
	  for (itr = L.begin(); itr != L.end(); itr++)
	  {
	  	i = *itr;
	  	
	  	double d = abs(A[i][j]);
	  	
	  	// dividing all elements of row i by absolute coeeficient of A[i][j]
	  	for(int k= 0; k<columns;k++){
	  		A[i][k] = A[i][k] / d;
	  	}
	  	b[i]=b[i]/d;
	  }
	  }
//	   printmatrix(A);
	   if( !U.empty()){
	   for (itr = U.begin(); itr != U.end(); itr++)
	  {
	  	i = *itr;
	  	double d = abs(A[i][j]);
	  	
	  	// dividing all elements of row i by absolute coeeficient of A[i][j]
	  	for(int k= 0; k<columns;k++){
	  		A[i][k] = A[i][k] / d;
	  		
	  	}
	  	b[i]=b[i]/d;
	  	
	    
	  }
	  }
//	   printmatrix(A);
 //write exit condition
	  if(columns == 1) //indicates last variable
	  {
		if(!L.empty()){
		vector<double> lowB ;
		 for (itr = L.begin(); itr != L.end(); itr++)
	 	 {
	  		i = *itr;
	  		lowB.push_back(b[i]);
	  	}
	  	
		lb = (int)(-1)*floor(get_max(lowB));
		}
		if(!U.empty()){
		vector<double> upB ;
		 for (itr = U.begin(); itr != U.end(); itr++)
	 	 {
	  		i = *itr;
	  		upB.push_back(b[i]);
	  	}
		
		ub=(int)floor(get_min(upB));
		}
		tempRowA.clear();
		for(int i = lb; i <= ub; i++){
		tempRowA.push_back(i);
		sol.push_back(tempRowA);
		tempRowA.clear();
		}
		
//		printmatrix(sol);
		return sol;
	}

		// create new matrix with combined equations after eliminating one variable
		int new_rows = L.size() * U.size() + Z.size();
		int new_columns = columns - 1;
		i=0;
		k=0;
		tempRowA.clear();
		for (itr = L.begin(); itr != L.end(); itr++){
			i = *itr;
			for(itr1 = U.begin(); itr1 != U.end(); itr1++){
				k= *itr1;
				// add row Ai + row Ak 
				// add  Bi + Bk
				for(int l =0;l<columns-1; l++){
				tempRowA.push_back(A[i][l]+A[k][l]);
				}
				newA.push_back(tempRowA);
				tempRowA.clear();
				newb.push_back(b[i]+b[k]);
			}
		}
//		printmatrix(newA);
		if(Z.size() > 0){
		for(itr = Z.begin(); itr != Z.end(); itr++ ){
		int z = *itr;
			for(int l =0;l<columns-1; l++){
         			tempRowA.push_back(A[z][l]);	
				}
				newA.push_back(tempRowA);
				newb.push_back(b[z]);
		}
		}
//		printmatrix(newA);
		vector<vector<double>> solutions ;
		solutions = FME(new_rows,new_columns,newA,newb);
		
		//indicates no solution to original set
		if(solutions.empty()){
		cout<<"No Solution"<<endl; 
		return solutions;
		}
		return solutions;
	
}


int main(int argc, char **argv){
 if (argc != 2)
    {
        cout << "Usage: " << argv[0] << " FILENAME\n";
        return 1;
    }
    ifstream input(argv[1]);
    if(!input){
        cerr << "error opening file " << argv[1] << "\n";
        return 1;
    }

    int m = 0; // #rows
    int n = 0; // #columns

    input >> m >> n;

    vector<vector<double>> A(m, vector<double>(n)); // get A
    for (int i = 0; i< m; i++){
        for (int j = 0; j< n; j++){
            input >> A[i][j];
        }
    }
    vector<double> b(m) ; // get b
    for (int i = 0; i< m; i++){
        input >> b[i];
    }
    cout<<".............................................................\n";
    printmatrix(A);
    cout<<".............................................................\n";
    cout<<".............................................................\n";
    printlist(b);
    cout<<".............................................................\n";

	vector<vector<double>> s = FME(m,n,A,b);
	// print out vector
	if(s.size()>0)
		cout<<"Has solution\n";
	else
		cout<<"Does not have solution\n";
	return 0;
}



