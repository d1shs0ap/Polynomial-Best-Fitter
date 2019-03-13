#include <bits/stdc++.h>

// The construction of the best-fitting line uses the theorem A^TAx=A^ty
// Written by d1shs0ap from scratch
// References: linear algebra book by Anton and Rorres



using namespace std;

//transpose
vector<vector<double> > transpose(vector<vector<double> > matrix){
    
    vector<vector<double> > matrixT(matrix[0].size());
    
    for(int j = 0; j < matrix[0].size(); j++){
        for(int i = 0; i < matrix.size(); i++){
            matrixT[j].push_back(matrix[i][j]);
        }
    } 
    
    return matrixT;
}

//matrix multiplication (naive)


double dotProduct(vector<double> v1, vector<double> v2){
    
    double sum = 0;
    
    for(int i = 0; i < v1.size(); i++){
        sum += v1[i]*v2[i];
    }
    
    return sum;
}

vector<vector<double> > multiply(vector<vector<double> > matrix1, vector<vector<double> > matrix2){
    
    vector<vector<double> > product(matrix1.size());
    
    for(int i = 0; i < matrix1.size(); i++){
        for(int j = 0; j < matrix2[0].size(); j++){
            
            vector<double> row;
            vector<double> column;
            
            row = matrix1[i];
            
            for (int k = 0; k < matrix2.size(); k++){ //iterate through rows
                column.push_back(matrix2[k][j]);
            }
            
            double element = dotProduct(row, column);//take the dot M1 row and M2 column to find element
            
            product[i].push_back(element);
        
        }
    }
    
    return product;
}

//Gaussian Elimination method (of square matrix)

//the 1st argument is Augmented matrix, 2nd is # of rows function returns the answer
vector<double> gaussianElim(vector<vector<double> > Aug, int M){
    //iterate through all columns, forward elimination
    for(int j=0; j<M-1; j++){
        /*switch rows to the ones w the largest column value to ensure the
        entry is non-zero and lower errors by avoiding to divide
        small numbers by large numbers (which may cause a rounding error)*/
        
        int l = j;//keeps track of row w largest column value
        
        /*only look for rows below cuz rows above are 
        already in REF*/
        
        for(int i=j+1; i<M; i++){ //finds row w/ largest column value
            if(abs(Aug[i][j]) > abs(Aug[l][j])){//don't forget to take absolute value
                l = i;
            }
        }
        
        //switches the rows
        
        for(int k = 0; k<=M; k++){//goes through elements of a row
            
            double t;//place holder to perform [l][k]=[j][k], [j][k]=[l][k]
            t = Aug[l][k];
            Aug[l][k]=Aug[j][k];
            Aug[j][k]=t;
            
        }
        
        //now REF by subtracting a constant times the jth row
        
        for(int i = j+1; i<M; i++){//iterates through rows
            
            double c = Aug[i][j]/Aug[j][j];//constant required to set i,j zero
            
            for(int k = j; k<= M; k++){ //goes through elements of a row
               Aug[i][k] -= c*Aug[j][k]; //what happens when k=j?
            }
            
        }
        
    }
    
    double sol[M]; //here is the solution
    
    //back substitution(we now have an upper triangular matrix)
    for(int i=M-1; i>=0; i--){//iterates through rows starting from the last(why?)
        
        double t = 0; //collects the sum of the rest of LHS = a_(i+1)(k+1)x_(k+1)+...
        
        for(int k = i+1; k < M; k++){//only interested in the upper triangle, i.e. k>i
            t += Aug[i][k]*sol[k]; //a_(i+1)(k+1)x_(k+1)
        }
        
        
        sol[i] = (Aug[i][M]-t)/Aug[i][i];/*subtract from the rightmost column(constant column),
                                 then divide by the entry to get x_i (why?)*/
        
    }
    /* Now change the solution which was in array form to vector because:
     * 1. Array allows me to assign values directly instead of push_back
     * 2. An array cannot be returned in a function */
    vector<double> solution;
    
    for(int i=0; i<M; i++){
        solution.push_back(sol[i]);
    }
    
    return solution;
}

//power function
double power(double base, int exponent){
    double product = 1;
    for(int i = 0; i<exponent; i++){
        product *= base;
    }
    return product;
}


int main() {
    //first ask for number of points
    printf("How many points?");
    int N;
    scanf("%d", &N);
    
    //what are the points?
    
    vector<int> xCoordinates;
    vector<int> yCoordinates;
    
    for (int i = 1; i<=N; i++){
        double x;
        double y;
        
        scanf("%lf %lf", &x, &y);
        
        xCoordinates.push_back(x);
        yCoordinates.push_back(y);
        
    }
    
    //ask for degree of polynomial
    printf("Enter degree of polynomial: ");
    int deg;
    scanf("%d", &deg);
    
    /*
    for(auto x: xCoordinates){
        cout << x <<"\n";
    }*/
    
    //create A
    vector<vector<double> > A (xCoordinates.size());
    for(int i = 0; i < xCoordinates.size(); i++){ //row
        for(int j = 0; j<=deg; j++){ //columns
            A[i].push_back(power(xCoordinates[i], j));
            
        }
    }
    
    
    
    //Find A^T
    vector<vector<double> > AT = transpose(A);
    
    
    //Find A^TA
    vector<vector<double> > matrixLHS = multiply(AT,A);
    
    /*for(auto row: matrixLHS){
        for(auto element: row){
            cout << element << " ";
        }
        cout << "\n";
    }*/
    
    
    //create y
    vector<vector<double> > Y(yCoordinates.size());
    for(int i = 0; i < yCoordinates.size(); i++){
        Y[i].push_back(yCoordinates[i]);
    }
    
    //Find A^Ty
    vector<vector<double> > matrixRHS = multiply(AT, Y);
    
        
    
    //Create augmented matrix
    vector<vector<double> > Aug = matrixLHS;
    for(int i = 0; i <= deg; i++){ //rows
        Aug[i].push_back(matrixRHS[i][0]);
    }
    
    /*for(auto row: Aug){
        for(auto element: row){
            cout << element << " ";
        }
        cout << "\n";
    }*/
    
    
    //Solve A^TAx=A^Ty via Gaussian Elimination
    vector<double> X = gaussianElim(Aug,deg+1);
    
    for(auto x: X){
        cout << x << " ";
    }
    
}

