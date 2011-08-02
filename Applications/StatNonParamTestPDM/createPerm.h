#ifndef CREAT_PERM_H__MST
#define CREAT_PERM_H__MST

#define GROUP_A_LABEL   0
#define GROUP_B_LABEL   1

typedef struct
  {
  int randNum;
  int index;
  } PermElement;

void generatePerm(int length, int * genPerm);

// generate a permutation of numbers from 0 to length-1
void generatePermGroup(int * groupID, int lengthGroupA, int lengthGroupB, int * genGroupID);

// generate a permutation of group assignments

int smallerPermElem(PermElement * elem1, PermElement * elem2);

#endif
