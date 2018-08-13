
#include <math.h>
#include <ctype.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>

using namespace std;

const int MaxBaseSize=1000;
const bool PRINTROWS=true;

int NumInputs; 
int DepthLimit;
int Threshold;
int MiniDistance;
int NumTargets;
int ProgramSize;
long long int Target[MaxBaseSize];
long long int NewIndex[MaxBaseSize];
int Dist[MaxBaseSize]; //distance from current base to Target[i]
int NDist[MaxBaseSize]; //what Dist would be if NewBase was added
long long int Base[MaxBaseSize];
int BaseSize;
int TargetsFound;
char Result[MaxBaseSize][50];
int  Res;
char *flag;
int Depth[MaxBaseSize];
int MaxDepth;
int Cou;

int bitcount(long long int x);
int ReturnLoop();
void ChangeMatrix();
void InitBase();
void ReadTargetMatrix();
bool is_target(long long int x);
bool is_base(long long int x);
int NewDistance(int u); //calculates the distance from the base to Target[u]
int TotalDistance(); //returns the sum of distances to targets
bool reachable(long long int T, int K, int S);
bool EasyMove(); //if any two bases add up to a target, pick them
void PickNewBaseElement();
void binprint(long long int x); //outputs last NumInputs bits of x 

void PrintExpression(int No);
void PrintBase(void);
void PrintResult();
int Max(int a, int b);

ifstream TheMatrix;
//ofstream out_file;

int
main(int argc, char *argv[])
{
  int NumMatrices; 
  int sum;
  int order[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  int order1_2[16] = {0, 1, 6, 7, 2, 3, 4, 5, 8, 9, 14, 15, 10, 11, 12, 13};
  int order1_3[16] = {1, 0, 5, 4, 6, 7, 2, 3, 9, 8, 13, 12, 14, 15, 10, 11};
  int order2_2[16] = {1, 0, 6, 7, 3, 2, 4, 5, 9, 8, 14, 15, 11, 10, 12, 13};
  int order2_3[16] = {1, 0, 4, 5, 7, 6, 2, 3, 9, 8, 12, 13, 15, 14, 10, 11};
  int U[16], B[16], i;

  clock_t t1=clock();
  //TheMatrix.open("matrices/matrices.txt");
  TheMatrix.open("MDS_84XOR.txt");
  //TheMatrix.open("tlt");
  //TheMatrix.open("baaaaad");
  //TheMatrix.open("matrices/b");
  //TheMatrix.open("matrices.10.15");
  TheMatrix >> NumMatrices;
  TheMatrix >> Threshold;
  TheMatrix >> DepthLimit;
  for (i = 0; i < NumMatrices; i++)
  {
    ReadTargetMatrix(); 
    TargetsFound = 0;
    InitBase();
    ProgramSize = 0; 
    int counter = 0; 
    while (TargetsFound < NumTargets-Cou) 
    {
	  counter++;
	  if(ProgramSize > Threshold)
	  {
		break;
	  }
	  if (!EasyMove()) PickNewBaseElement();
    }
    if(ProgramSize <= Threshold)
    	PrintResult();
    else
    	printf("Out of threshold\n\n");
//    if(i % 3 == 2)
//    	cout << endl;

//    if(i % 2 == 0)
//    	sum = ProgramSize;
//    else
//    	if(sum + ProgramSize > 0)
//    	{
//    		cout << (i-1)/2 << endl;
//    		cout << sum << ' ' << ProgramSize << endl << endl;
//		}
	
//	if(i % 2 == 0)
//		U[i/2] = ProgramSize;
//	else
//		B[i/2] = ProgramSize;    	
  }
//  for(i = 0; i < 16; ++i)
//  {
//  	int k = order2_3[i];
//  	printf("%d\n%d\t%d\n\n", k, U[k], B[k]);
//  }
  	
  clock_t t2=clock();
  cout << (t2-t1)/(double)CLOCKS_PER_SEC;
}//main

int ReturnLoop()
{
    MiniDistance = BaseSize*NumTargets;
    //clock_t t1=clock();
    ChangeMatrix(); 
    //printf("%d \n", TargetsFound);
    InitBase();
    //printf("%d \n", TargetsFound);
    //ProgramSize = 0; 
    int counter = 0; 
    //printf("%d %d \n", TargetsFound, NumTargets);
    while (TargetsFound < NumTargets-Cou) 
    {
	  counter++;
	  //printf("%d \n", 12345);
	  if(ProgramSize > Threshold)
	  {
		break;
	  }
	  //printf("%d \n", 12345);
	  if (!EasyMove()) PickNewBaseElement();
	  //printf("%d %d \n", TargetsFound, NumTargets-Cou);
      
    }
    if (ProgramSize > Threshold)
    	printf("Out of threshold\n\n");
    	
  
  //clock_t t2=clock();
  //cout << (t2-t1)/(double)CLOCKS_PER_SEC;
}

void InitBase()
{
	
	Res = ProgramSize;
	TargetsFound = 0;
	Base[0] = 1;
	Depth[0] = 0;
	MaxDepth = 0;
	Cou = 0;
	for (int i = 1; i < NumInputs; i++) 
	{
		Base[i] = 2*Base[i-1];
		Depth[i] = 0;
	}
	BaseSize = NumInputs; //initial base is just the xi's, depth is 0
	for (int i = 0; i < NumTargets; i++) {
	
	    if (Target[i] == 0) Cou += 1; continue;
	    if (Dist[i] == 0) 
		{
			TargetsFound++;
			//printf("%d \n", i);
			//print the expression of output and input
			for(int j = 0; j < NumInputs; ++j)
				if(Base[j] == Target[i])
				{
					sprintf(Result[Res], "y%d = x%d  *  (0)\n", i, j);
					++Res;
					break;
				}
		}
	}
}

int TotalDistance() //returns the sum of distances to targets
{
  int D = 0;
  int t;
  for (int i = 0; i < NumTargets; i++) 
  {
    t = NewDistance(i); 
    NDist[i] = t; 
    D = D + t;
  }
  return D;
}

long long int NewBase; //global variable containing a candidate new base

bool EasyMove()
{
  int t;
  bool foundone = false;
  
  //see if anything in the distance vector is 1
  for(int i = 0; i < NumTargets; i++) 
    if (Dist[i] == 1) 
    {
      foundone = true;
      t = i;
      //printf("y%d \n", t);
      MiniDistance = 0;
      break;
    }
  if (!foundone) return false;
  //update Dist array
  NewBase = Target[t]; 
  for (int u = 0; u < NumTargets; u++) 
  {
    Dist[u] = NewDistance(u); 
	MiniDistance += Dist[u]; 
  }
  //update Base with NewBase 
  Base[BaseSize] = NewBase; 
  BaseSize++; 
  ProgramSize++;
  //print the expression
  PrintExpression(t);
  TargetsFound++;
  return true;
} //EasyMove()

void PrintResult()
{
	cout << ProgramSize << endl;
	for(int i = 0; i < Res; ++i)
		printf("%s", Result[i]);
	cout << "Depth is " << MaxDepth << endl << endl;
}

int Max(int a, int b)
{
	if(Depth[a] > Depth[b])
		return Depth[a];
	else
		return Depth[b];
}

/* print the expression*/
void PrintExpression(int No)
{
	int i, j;
	
	for(i = 0; i < BaseSize; ++i)
		for(j = i + 1; j < BaseSize; ++j)
			if((Base[i] ^ Base[j]) == Base[BaseSize-1])
			{
				Depth[BaseSize-1] = Max(i, j) + 1;
				if(Depth[BaseSize-1] > MaxDepth)
					MaxDepth = Depth[BaseSize-1];
				flag = Result[Res];
				flag += sprintf(flag, "t%d = ", ProgramSize);
				if(i < NumInputs)
					flag += sprintf(flag, "x%d + ", i);
				else
					flag += sprintf(flag, "t%d + ", i - NumInputs + 1);
				if(j < NumInputs)
					flag += sprintf(flag, "x%d *  y%d  (%d) %#x\n", j, No, Depth[BaseSize-1], Base[BaseSize-1]);
				else
					flag += sprintf(flag, "t%d *  y%d  (%d) %#x\n", j - NumInputs + 1, No, Depth[BaseSize-1], Base[BaseSize-1]);
				++Res;
				return;
			}
}

/* PickNewBaseElement is only called when there are no 1's in Dist[]*/
void PickNewBaseElement()
{
	int MinDistance;
	long long int TheBest;
	int ThisDist;
	int ThisNorm, OldNorm;
	int besti,bestj, d;
	bool easytarget;
	int BestDist[MaxBaseSize];
	
	MinDistance = BaseSize*NumTargets; //i.e. something big
	OldNorm = 0; //i.e. something small
	//try all pairs of bases
	for (int i = 0; i < BaseSize - 1; i++)
	{
		if (Depth[i] + 1 >= DepthLimit) continue;
		for (int j = i+1; j < BaseSize; j++)
		{
			if (Depth[j] + 1 >= DepthLimit) continue;
			NewBase = Base[i] ^ Base[j];
			//sanity check
			if (NewBase == 0) { cout << "a base is 0, should't happen " << endl; exit(0); }
			//if NewBase is not new continue
			if (is_base(NewBase)) continue;
			//if NewBase is target then choose it
			easytarget = false;
			if (is_target(NewBase))
			{
				cout << "shouldn't find an easy target here " << endl; 
				exit(0);
				easytarget = true;
				besti = i;
				bestj = j;
				TheBest = NewBase;
				break;
			}
			ThisDist = TotalDistance(); //this also calculates NDist[]
			if (ThisDist <= MinDistance)
			{
				//calculate Norm
				ThisNorm = 0; 
				for (int k = 0; k < NumTargets; k++)
				{
					d = NDist[k];
					ThisNorm = ThisNorm + d*d;
				}
				//resolve tie in favor of largest norm
				if ((ThisDist < MinDistance) || (ThisNorm > OldNorm) )
				{
					besti = i;
					bestj = j;
					TheBest = NewBase;
					for (int uu = 0; uu < NumTargets; uu++) BestDist[uu] = NDist[uu]; 
					MinDistance = ThisDist;
					OldNorm = ThisNorm;
				}
			}   
		}
	if (easytarget) break;
	}
	if (MinDistance == MiniDistance) ReturnLoop();
	else {
	MiniDistance = MinDistance;
	//update Dist array
	NewBase = TheBest; 
	for (int i = 0; i < NumTargets; i++) Dist[i] = BestDist[i];
	//printf ("%d %d %d \n", besti, bestj, MiniDistance);
	//std::cout<<"Press [ENTER] to continue\n";
    //std::cin.get();
	//update Base with TheBest
	Base[BaseSize] = TheBest;
	Depth[BaseSize] = Max(besti, bestj) + 1; 
	if(Depth[BaseSize] > MaxDepth)
		MaxDepth = Depth[BaseSize];
	BaseSize++;
	//output linear program
	ProgramSize++;
	//print the expression
	flag = Result[Res];
	flag += sprintf(flag, "t%d = ", ProgramSize);
	if(besti < NumInputs)
		flag += sprintf(flag, "x%d + ", besti);
	else
		flag += sprintf(flag, "t%d + ", besti - NumInputs + 1);
	if(bestj < NumInputs)
		flag += sprintf(flag, "x%d  (%d) %#x\n", bestj, Depth[BaseSize-1], Base[BaseSize-1]);
	else
		flag += sprintf(flag, "t%d  (%d) %#x\n", bestj - NumInputs + 1, Depth[BaseSize-1], Base[BaseSize-1]);
	++Res;
	if (is_target(TheBest)) TargetsFound++; //this shouldn't happen
    }
} //PickNewBaseElement()

void binprint(long long int x) //outputs last NumInputs bits of x 
{
  long long int t = x;
  for (int i = 0; i < NumInputs; i++)
  {
    if (t%2) cout << "1 "; else cout << "0 ";
    t = t/2;
  }
} //binprint  

void PrintBase()
{
	int i;
	
	for(i = 0; i < BaseSize; ++i)
	{
		binprint(Base[i]);
		printf("\n");
	}
}

void ReadTargetMatrix()
{
  TheMatrix >> NumTargets;
  TheMatrix >> NumInputs;
  //check that NumInputs is < wordsize
  if (NumInputs >= 8*sizeof(long long int)) 
  {
    cout << "too many inputs" << endl;
    exit(0);
  }

  int bit;
  for (int i = 0; i < NumTargets; i++)
  //read row i
  {
    long long int PowerOfTwo  = 1;
    Target[i] = 0;
    Dist[i] = -1; //initial distance from Target[i] is Hamming weight - 1
    TheMatrix.ignore(1, '[');
    for (int j = 0; j < NumInputs; j++) 
    {
      TheMatrix >> bit;
      if (bit) 
      {
        Dist[i]++; 
		Target[i] = Target[i] + PowerOfTwo;
      }
      PowerOfTwo = PowerOfTwo * 2;
    }
	TheMatrix.get();
  }
} //ReadTargetMatrix()

void ChangeMatrix()
{
  
  for (int i = 0; i < NumTargets; i++)
  {
  	
	if (Dist[i] == 0) {
	    
	    Dist[i] = -1;
	    Target[i] = 0;
    }
	else Dist[i] = bitcount(Target[i]) - 1;
  }
  
}

bool is_target(long long int x)
{
  for (int i = 0; i < NumTargets; i++)
    if (x == Target[i]) return true;
  return false;
} //is_target

bool is_base(long long int x)
{
  //sanity check, shouldn't ask if 0 is base
  if (x == 0) {cout << "asking if 0 is in Base " <<endl ; exit(0); }
  for (int i = 0; i < BaseSize; i++) if (x == Base[i]) return true;
  return false;
} //is_base

// Distance is 1 less than the number of elements
// in the base that I need to add in order to get Target[u].
// The next function calculates the distance from the base,
// augmented by NewBase, to Target[u]. Uses the following observations:
// Adding to the base can only decrease distance. 
// Also, since NewBase is the sum of two old base 
// elements, the distance from the augmented base 
// to Target[u] can decrease at most by 1. If the
// the distance decreases, then NewBase must be one
// of the summands.
  
int NewDistance(int u) 
{
  //if Target[u] is in augmented base return 0;
  if (Target[u] == 0) return 0;
  else 
  {
    if (is_base(Target[u]) || (NewBase == Target[u])) return 0;
  
  // Try all combinations of Dist[u]-1 base elements until one sums 
  // to Target[u] + NewBase. If this is true, then Target[u] is the
  // sum of Dist[u] elements in the augmented base, and therefore
  // the distance decreases by 1.
  
    if (reachable(Target[u] ^ NewBase,Dist[u]-1,0)) return (Dist[u]-1);
    else return Dist[u]; //keep old distance 
  }
} //NewDistance(int u) 


//return true if T is the sum of K elements among Base[S..BaseSize-1]
bool reachable(long long int T, int K, int S)
{
    if ((BaseSize-S) < K) return false; //not enough base elements
    
    if (K==0) return false; //this is probably not reached
    if (K==1) 
    {
      for (int i=S; i < BaseSize; i++) if (T == Base[i]) return true;
      return false;
    } 
    
    //consider those sums containing Base[S]
    if (reachable(T^Base[S], K-1, S+1)) return true;
    //consider those sums not containing Base[S]
    if (reachable(T, K, S+1)) return true;
    //not found
    return false;
} // reachable(long long int T, int K, int S)

int bitcount(long long int x)
{
	int C = 0;
	while (x) {
		if (x & 0x1) C++;
		x >>= 1;
	}
	return C;
}
