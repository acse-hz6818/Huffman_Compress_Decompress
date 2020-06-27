#pragma once
#include<iostream>
#include<fstream>
using namespace std;

int leafArrIndex = 0; 
int strIndex = 0;
int nodeIndex = 0;

class TreeNode
{
public:
	TreeNode(){
		parent = lchild = rchild = -1;
	}
	TreeNode(const char &data, const int &wt);
	virtual ~TreeNode(){}
public:
	char value;
	int weight;
	int parent, lchild, rchild;

};

TreeNode::TreeNode(const char &data, const int &wt){
	value = data;
	weight = wt;
	parent = lchild= rchild = -1;
}


class HuffmanCode{
public:
	HuffmanCode() : length(20){
		ptr = new char[length];
	}
	virtual ~HuffmanCode(){
		if (ptr)
		{
			delete[] ptr;
			ptr = NULL;
		}
	}
public:
	char *ptr;
	const int length;
};

class HuffmanTree{
public:
	HuffmanTree(){
		maxSize = 1000, 
		arrayTree = new TreeNode[maxSize];
		codeArray = NULL;
		currsize = 0;
	}
	virtual ~HuffmanTree() {
		if (codeArray)
		{
			delete[]codeArray;
			codeArray = NULL;
		}
		if (arrayTree)
		{
			delete[]arrayTree;
			arrayTree = NULL;
		}
	}

	void Process(const char*, const char*);
	void ReBuildHuffmanTree(const char*, const int*, int);
	void HuffmanCoding();
	void ReBuildArray(char **, char *, int);

public:
	int currsize;
	TreeNode *arrayTree;
	HuffmanCode *codeArray;

private:
	int maxSize;
	int FindPosition(const char &) const; 
	void Insert(const char&, const int&);
	void BuildHuffmanTree();
	void WwriteCodeToFile(const char*);
	void WriteTreeStructureToFile(const int &, ofstream&fw);
	void WriteLeafToFile(const int &, ofstream&fw);
};

