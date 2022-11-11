# CSE521 Assignment 1 : Fourier-Motzkin Elimination
 
Implemented Fourier Motzkin Elimination (FME) in three steps, as follows: </br>
## (a) Implemented in BaselineFME.cpp. Takes as input a system of inequalities, Ax ≤ b, where A is a constant matrix of mxn (m: number of inequalities and n: number of unknowns), b is an m-entry constant vector and x is an n-entry vector of unknowns. Implementation will output whether the system of inequalities have any solution or not.  
### Input File format
```
m n
a11 a12 a13 .. a1n 
a21 a22 a23 .. a2n
..
..
am1 am2 am3 .. amn
b1
b2
..
bm
```
where m is number of rows and n is number of unknowns

### Sample 
For Inequations such as this
```
x1 - 4x2 <=2
x1 + 5x2 <=7
-x1 <= -3
```
Input file should be
```
3 2
1 -4
1 5
-1 0
2
7
-3
```
### To run

```
$ g++ BaselineFME.cpp
$ ./a.out input.txt

```
### Sample Output for (1)
![image](https://user-images.githubusercontent.com/112839333/201240189-de8f2924-a167-4767-9585-d8965caea6de.png)


## (b) This is  an  integer  version  of  FME. 
 
## (b.1) In the first version, your program will declare that there is an integer solution if and only if all projections (reductions) are exact. In addition, your code should indicate whether  the projection  you employed at each step of the solution process was exact or inexact.   
### To run

```
$ g++ FME2aInt.cpp
$ ./a.out input.txt

```
### Sample Output for (2.b.1)
![image](https://user-images.githubusercontent.com/112839333/201255025-ca132e03-728a-4b03-9633-c17765c9be9e.png)

## (b.2) In the second version, your program is going to form “dark shadow” equations.  In both the integer versions, your implementation will print out a loop nest which, when executed, prints all the integer points in the solution space.  
### To run

```
$ g++ FME2bDark_shadows.cpp
$ ./a.out input.txt

```
### Sample Output for (2.b.2)
![image](https://user-images.githubusercontent.com/112839333/201260033-b2d5402d-e429-4865-8ce8-961c2b912cd5.png)


