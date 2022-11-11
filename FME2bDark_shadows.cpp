/*

ALgorithm FMInt
procedure projectint(A,b)
	inexact = False
	Repeat
		choose an unkown xj to project
		L = set of row indices i, where aij <0
		U = set of row indices i, where aij >0
		if L or U is empty, its an unconstrained variable.. exit
		for i that belongs to L union U do
			g=|gcd(ai1,ai2,...)|
			divide the whole row Ai by g
			let bi = floor bi/g
		end for
		for i belongs to L do
			for k belongs to U do
				add a new inequality akj* row Ai - aij * Row Ak <=  akj*bi - aij * bk + akj * aij + akj - aij -1 // darkshadow
				if |akj * aij| >1 then inexact = True
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
#include <list>  
#include <algorithm>  
#include <iterator>  
#include<bits/stdc++.h>
using namespace std;
typedef struct sols{
	vector<vector<double>> solutions;
	bool inexact;
	bool darkempty;
}sols;
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

int gcd(int a, int b)
{
    if (a == 0)
        return b;
    return gcd(b % a, a);
}

int findGCD(vector<double> array)
{
    int result = (int)array[0];
    for (int i = 1; i < array.size(); i++)
        result = gcd((int)array[i], (int)result);
  
    return result;
}
  

void printmatrix(vector<vector<double> > &A){
cout<<"Print"<<endl;
	for(size_t i=0; i<A.size(); i++){
		for(auto j: A[i]){
			cout<<j<<" ";
		}
		cout<<endl;
	}
	cout<<"Print end"<<endl;
	return;
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
void printset(set<int>  &s)
{ cout<<"Print set"<<endl;
    for (auto  &i: s) {
        cout << i << " ";
    }cout<<"Print set end"<<endl;
    return;
}
vector<vector<double>> subfMatrix(vector<vector<double>> &A,int index){
	vector<vector<double>> result;
	for(int i=0; i<A.size(); i++){
		result.push_back(vector<double> (A[i].begin()+index, A[i].begin()+index+1));
	}
	return result;
}
vector<vector<double>> sublMatrix(vector<vector<double>> &A){
	vector<vector<double>> result;
	for(int i=0; i<A.size(); i++){
		result.push_back(vector<double> (A[i].begin(), A[i].end()-1));
	}
	return result;
}
vector<vector<double> > multiply(vector<vector<double>>&A, vector<vector<double>>&B){
    vector<vector<double> >result(A.size(),vector<double>(B[0].size(),0));
    if(A[0].size()==B.size()){
        const int N=A[0].size();
        for(int i=0;i<A.size();i++){
            for(int j=0;j<B[0].size();j++){
                for(int k=0;k<A[0].size();k++)
                    result[i][j]+=A[i][k]*B[k][j];
            }
        }
    }
//    printmatrix(result);
    return result;
}
sols FME( int rows, int columns, vector<vector<double>> A,vector<double> b){
			
	// j is the first variable we choose to eliminate, let j = columns -1;
	int j = columns-1,i,k,lb =0,ub =0, gcd =1;
    sols  sol ;
	bool inexact = false;
	bool darkEmpty = false;
	set<int> L,U,Z;
	set<int>::iterator itr,itr1;
	vector<vector<double>> newA ;
	vector<double> newb,tempRowA;
	double akj,aij;
	for(int i=0;i<rows;i++){
		if(A[i][j] < 0){
			L.insert(i);
		}
		else if (A[i][j] > 0){
			U.insert(i);
		} else {
			Z.insert(i);
		}
	}
//	printset(L);
//	printset(U);
//	printset(Z);
	
	if( (L.empty() || U.empty()) ){//&& columns !=1){
		cout<<"unconstrained variable x"<<j+1<<endl<<"No closed iteration space"<<endl;
		return sol;
		
	}
	if( !L.empty() ){
	  for (itr = L.begin(); itr != L.end(); itr++)
	  {
	  	i = *itr;
	  	gcd = abs(findGCD(A[i]));
	  	
	  	// dividing all elements of row i by absolute gcd
	  	for(int k= 0; k<columns;k++){
	  		A[i][k] = A[i][k] / gcd;
	  		
	  	}
	  	b[i]=floor(b[i]/gcd);
	    
	  }
	  }
	   if( !U.empty()){
      for (itr = U.begin(); itr != U.end(); itr++)
	  {
	  	i = *itr;
	  	gcd = abs(findGCD(A[i]));
	  	// dividing all elements of row i by absolute coeeficient of A[i][j]
	  	for(int k= 0; k<columns;k++){
	  		A[i][k] = A[i][k] / gcd;
	  	}
	  	b[i]=floor(b[i]/gcd);
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
		lb = (int)(-1)*floor(get_min(lowB));
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
		sol.solutions.push_back(tempRowA);
		tempRowA.clear();
		}
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
				akj = A[k][j];
				aij = A[i][j];
				// add row Ai + row Ak 
				// add  Bi + Bk
				for(int l =0;l<columns-1; l++){
				tempRowA.push_back((akj*A[i][l])-(aij*A[k][l]));
				}
				newA.push_back(tempRowA);
				tempRowA.clear();
				//computing dark shadow
				newb.push_back((akj * b[i])-( aij * b[k]) + (akj * aij) + akj - aij -1);
				//if |akj * aij| >1 then inexact = True
				if (abs(akj * aij) > 1){
				inexact =true;
				sol.inexact =true;
				cout<<"Inexact Projection"<<endl;
				} 
				else 
				cout<<"Exact Projection"<<endl;
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
		sols ss= FME(new_rows,new_columns,newA,newb);
		solutions = ss.solutions;
		ss.inexact = sol.inexact;
		//indicates no solution to original set
		if(solutions.empty()){cout<<"solutions empty"; return sol;}
		// process solutions
		// Calculate new bounds on current variable
		int nsols = solutions.size();
		int vsolved = solutions[0].size();
		
		vector<double> temp ;
//		printmatrix(solutions);
		vector< vector<double> > matX;
		for(int jc =0;jc<solutions[0].size();jc++){//1
			for(int ir=0;ir< solutions.size() ; ir++){//4
				temp.push_back(solutions[ir][jc]);
			}
			matX.push_back(temp);
		}
//		printmatrix(matX);
		vector<vector<double>> matAX;
//		printmatrix(A);
		vector<vector<double>> subA = subfMatrix(A,vsolved);// strip only first column
//		printmatrix(subA);
		vector<vector<double>> subMA = sublMatrix(A);//strip last coumn
		matAX=multiply(subMA,matX);
//		printmatrix(matAX);
		sols finalsols ;
		for(int col = 0;col < nsols;col++){
		 vector<double> LBList ;
		 vector<double> UBList ;
//       printlist(b);
	   for (itr = L.begin(); itr != L.end(); itr++)
	  	{
	  	i = *itr;
		 LBList.push_back((b[i] - matAX[i][col])/subA[i][0]);
		 }
	    for (itr = U.begin(); itr != U.end(); itr++)
	  	{
	  	i = *itr;
		 UBList.push_back((b[i] - matAX[i][col])/subA[i][0]);
		 }
//         printlist(LBList);
//         printlist(UBList);
		int lbb=0,ubb=0;
		if(!LBList.empty())
		lbb =  (int)((-1)*floor(get_min(LBList)));	
      	if(!UBList.empty())
        ubb =(int)floor(get_min(UBList));
//	        printmatrix(solutions);
        
	for(int i = lbb;i <= ubb;i++){
		vector<double> nv(solutions[col]);
		nv.push_back(i); //inser at end
		finalsols.solutions.push_back(nv);
	}
	}
	finalsols.inexact = ss.inexact ;
	return finalsols;
	
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

sols s = FME(m,n,A,b);
cout<<"Done computation"<<endl;
 // print out vector
 printmatrix(s.solutions);
 

return 0;

}



