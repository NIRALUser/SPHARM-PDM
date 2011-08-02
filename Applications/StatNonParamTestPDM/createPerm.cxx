#include <algorithm>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <ctime>

#include "createPerm.h"

using namespace std;

static int first = 0;

/* Inquire with Brad whether he has already written a permutation class, otherwise write one
template <VariableType> permuteVariable( VariableType originalVariables, VariableType permutedVariable)
// generate a permutation of the supplied vector of variables
// VariableType has to be a vector/matrix type that allows assignment and Size() defined
{
  int numSubjects = originalVariables.Size();
  if (!first) {
    first = 1; srand(time(NULL));
    //cout << "generatePermGroup called" << endl;
  }
  int * newPerm = new int [numSubjects];
  generatePerm(numSubjects, newPerm);
  for (int i=0; i<numSubjects; i++) {
    permutedVariable[i] = originalVariables[newPerm[i]];
  }

  delete newPerm;
}
*/

void generatePermGroup(int * groupID, int lengthGroupA, int lengthGroupB,
                       int * genGroupID)
// generate a permutation of group assignments
{
  int numSubjects = lengthGroupA + lengthGroupB;

  if( !first )
    {
    first = 1; srand(time(NULL) );
    // cout << "generatePermGroup called" << endl;
    }
  int * newPerm = new int[numSubjects];
  generatePerm(numSubjects, newPerm);
  for( int i = 0; i < numSubjects; i++ )
    {
    genGroupID[i] = groupID[newPerm[i]];
    }

  delete newPerm;

}

void OLDgeneratePermGroup(int * groupID, int lengthGroupA, int lengthGroupB,
                          int * genGroupID)
// generate a permutation of group assignments
{
  int numSubjects = lengthGroupA + lengthGroupB;
  int numMaxSwaps = lengthGroupA;

  if( lengthGroupB < numMaxSwaps )
    {
    numMaxSwaps = lengthGroupB;
    }

  // Replace drand48() - Matthieu March,02 2004
  int numSwaps = (int) ( (float) (double(rand() ) / RAND_MAX) * (numMaxSwaps / 2 - 1) + numMaxSwaps / 2 );

  int * genPermA = new int[lengthGroupA];
  int * genPermB = new int[lengthGroupB];

  generatePerm(lengthGroupA, genPermA);
  generatePerm(lengthGroupB, genPermB);
  for( int i = 0; i < numSubjects; i++ )
    {
    genGroupID[i] = groupID[i];
    }
  // the first numSwaps encounters of Group A (scrambled according to genPermA) are
  // swaped with the those of Group B (scrambled according to genParamB)
  for( int swap = 0; swap < numSwaps; swap++ )
    {
    int indexA = 0, indexB = 0;
    int numCur = 0;
    int i;
    for( i = 0; i < numSubjects; i++ )
      {
      if( groupID[i] == GROUP_A_LABEL )
        {
        if( numCur == genPermA[swap] )
          {
          indexA = i;
          }
        numCur++;
        }
      }
    for( i = 0; i < numSubjects; i++ )
      {
      if( groupID[i] == GROUP_B_LABEL )
        {
        if( numCur == genPermB[swap] )
          {
          indexB = i;
          }
        numCur++;
        }
      }
    genGroupID[indexA] =  GROUP_B_LABEL;
    genGroupID[indexB] =  GROUP_A_LABEL;
    }
  delete genPermA;
  delete genPermB;

}

void generatePerm(int length, int * genPerm)
{
  if( !first )
    {
    first = 1; srand(time(NULL) );
    // cout << "generatePerm called" << endl;
    }
  PermElement * newPerm = new PermElement[length];
  int           cnt;
  for( cnt = 0; cnt < length; cnt++ )
    {
    newPerm[cnt].randNum = rand();
    newPerm[cnt].index = cnt;
    }
  qsort(newPerm, length, sizeof(PermElement),
        (int (*)(const void *, const void *) )smallerPermElem);
  for( cnt = 0; cnt < length; cnt++ )
    {
    genPerm[cnt] = newPerm[cnt].index;
    }
  delete newPerm;
}

int smallerPermElem(PermElement * elem1, PermElement * elem2)
{
  if( elem1->randNum > elem2->randNum )
    {
    return 1;
    }
  else if( elem1->randNum < elem2->randNum )
    {
    return -1;
    }
  else
    {
    return 0;
    }
}
