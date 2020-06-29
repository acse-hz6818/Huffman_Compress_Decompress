#include"HuffmanCompress_Un.h"



/*realization of Process function  ÷*/
void HuffmanTree::Process(const char *inFilename, const char *outFilename)
{
	ifstream inFile(inFilename, ios::in | ios::binary);       //Read in file to be compressed,  read binary binary file
	if (!inFile) {
		cerr << " Can not open£°" << endl;
		exit(1);
	}
	char ch = inFile.get();   //read one character

	while (!inFile.eof()) {   // read in characters from file until the ending, and count the number of characters
		int pos = FindPosition(ch);
		if (pos != -1)
			arrayTree[pos].weight++;
		else
			Insert(ch, 1);        //insert a new one, with weight 1

		ch = inFile.get();
	}
	BuildHuffmanTree();                   //build huffman tree
	HuffmanCoding();               // encoding for counted characters
	WriteCodeToFile(outFilename); // save code table into file
	inFile.close();
}


/*find position of ch  in the arrayTree[]*/
int HuffmanTree::FindPosition(const char &ch) const {
	for (int i = 0; i < currsize; i++)
		if (arrayTree[i].value == ch)
			return i;
	return -1;
}


/*insert TreeNode*/
void HuffmanTree::Insert(const char &data, const int &wt) {
	if (2 * currsize - 1 > maxSize) return;
	arrayTree[currsize].value = data;
	arrayTree[currsize].weight = wt;
	currsize++;
}


/*for compressing, build huffmanTree*/             //okokokokokok
void HuffmanTree::BuildHuffmanTree() {
	int i = currsize;
	int wt1, wt2;
	int s1 = 0, s2 = 0;
	while (i < 2 * currsize - 1) {
		wt1 = wt2 = numeric_limits<int>::max();     //returns the maximum value of double types allowed by the compiler
		for (int k = 0; k < i; k++) {
			if (arrayTree[k].weight < wt1 && arrayTree[k].parent == -1) {
				wt2 = wt1;
				s2 = s1;
				wt1 = arrayTree[k].weight;                      //set wt1 as minimal weight
				s1 = k;
			}
			else if (arrayTree[k].weight < wt2 && arrayTree[k].parent == -1) {
				wt2 = arrayTree[k].weight;                      //set wt2 as second minimal weight
				s2 = k;
			}
		}
		arrayTree[i].weight = arrayTree[s1].weight + arrayTree[s2].weight;
		arrayTree[i].lchild = s1;
		arrayTree[i].rchild = s2;
		arrayTree[s1].parent = arrayTree[s2].parent = i;
		i++;
	}
}

/*get exact code for each node*/    
void HuffmanTree::HuffmanCoding()
{
	char tempCodeArray[100] = { 0 };                   //store temporary code (initially inverted code)
	if (codeArray)
	{
		delete[]codeArray;
		codeArray = NULL;
	}
	codeArray = new HuffmanCode[2 * currsize - 1];   // store codes for each character
	if (!codeArray)
		return;

	int Len = 2 * currsize - 1;                           //total number of node
	int p, tmp_index;
	for (int i = 0; i < Len; i++) {
		tmp_index = i;
		p = arrayTree[i].parent;
		int n = 0;
		while (p != -1 && p < Len) {
			if (arrayTree[p].lchild == tmp_index)
				tempCodeArray[n++] = '0';         //add 0 if it's 0
			else if (arrayTree[p].rchild == tmp_index)
				tempCodeArray[n++] = '1';         //add 1 if it's 1

			tmp_index = p;
			p = arrayTree[tmp_index].parent;
		}
		for (int t = 0; t < n; t++) {            // Make codearray store codes with correct order
			codeArray[i].ptr[t] = tempCodeArray[n - 1 - t];

		}
		codeArray[i].ptr[n] = '\0';
	}
}


/*for compressing,  write HuffmanTree into file*/
void HuffmanTree::WriteCodeToFile(const char *outFilename) {

	ofstream fw(outFilename, ios::out | ios::binary);
	if (fw)
	{
		fw << currsize * 2 - 1 << " ";                   // write the total number of nodes
		WriteTreeStructureToFile(currsize * 2 - 2, fw);//start from root node, save rooNode in the end of arrayTree
		fw << " ";
		WriteLeafToFile(currsize * 2 - 2, fw);      //Traverse from the root node, write into ASCII code corresponding to leafNode character
		fw.close();
	}

}

/*for compressing, write HuffmanTree stducture into file*/
void HuffmanTree::WriteTreeStructureToFile(const int &num, ofstream &fw) {
	if (arrayTree[num].lchild == -1) {
		fw << "0";                         // leafnode is 0
		return;
	}
	else
		fw << "1";                         // internal node is 1

	WriteTreeStructureToFile(arrayTree[num].lchild, fw);
	WriteTreeStructureToFile(arrayTree[num].rchild, fw);
}


/*for compressing, write leafNode code into file*/
void HuffmanTree::WriteLeafToFile(const int &num, ofstream&fw)
{
	if (arrayTree[num].lchild == -1) {          // Only the left subtree is needed to be judged bcs it's full binary tree
		fw << (int)(arrayTree[num].value) << " "; //ASCII code corresponding to leafNode character
		return;
	}
	WriteLeafToFile(arrayTree[num].lchild, fw);
	WriteLeafToFile(arrayTree[num].rchild, fw);
}

/*compress file˛*/           
void CompressFile(const char *sourceFilename, const char *codeFilename, const char *geneFilename) {
	HuffmanTree tree;
	tree.Process(sourceFilename, codeFilename);

	ifstream inFile(sourceFilename, ios::in | ios::binary);
	ofstream outFile(geneFilename, ios::out | ios::binary);
	if (!inFile || !outFile) {
		cerr << "Open file fail!" << endl;
		exit(1);
	}

	int decimal = 0;                            //convert binary to decimal
	int i, j, k = 0;
	int bits = 0;                               //character encoding digits number
	int temp[10000];
	memset(temp, 0, sizeof(temp));


	while (!inFile.eof()) {
		decimal = 0;
		unsigned char one_byte = inFile.get();
		for (i = 0; i < tree.currsize; i++) {
			if (one_byte == tree.arrayTree[i].value) {
				bits += strlen(tree.codeArray[i].ptr);
				int len = strlen(tree.codeArray[i].ptr);  //record each codeArray[i] actual length
				if (bits < 8) {
					for (j = 0; j < len; j++)
						temp[k++] = tree.codeArray[i].ptr[j] - '0';
				}
				else if (bits >= 8) {
					for (j = 0; k < 8; j++)
						temp[k++] = tree.codeArray[i].ptr[j] - '0';
					for (; j < len; j++)
						temp[k++] = tree.codeArray[i].ptr[j] - '0';

					decimal += temp[0] * 128 + temp[1] * 64 + temp[2] * 32 + temp[3] * 16 + temp[4] * 8 + temp[5] * 4 + temp[6] * 2 + temp[7] * 1;

					for (j = 0; j < 8; j++)
						temp[j] = 0;                  //set 8-bit in temporary array to 0,  convenient for later data to use the array space

					for (j = 8; j < k; j++)
						temp[j - 8] = temp[j];        //将大于8位的部分向前移动, 用去除8位bit之后的数据填充新数组的前几位
					k = bits = j - 8;

					unsigned char c = decimal;   //discard the part less than 8 times, to achieve compression
					outFile.write((char*)&c, sizeof(c));
				}
			}
		}
	}
	if (bits) {             //after the end of the cycle, there may be  bit that is not enough 8-bit
		decimal += temp[0] * 128 + temp[1] * 64 + temp[2] * 32 + temp[3] * 16 + temp[4] * 8 + temp[5] * 4 + temp[6] * 2 + temp[7] * 1;
		unsigned char c = decimal;
		outFile.write((char*)&c, sizeof(c));
	}
	inFile.close();
	outFile.close();
	cout << "Compression Successful!" << endl;
}



/*for decompressing, writes information into array, consistent with the compression format*/
void HuffmanTree::ReBuildArray(char **HC, char *value, int Num) {
	int j = 0;
	for (int i = 0; i < Num; i++) {
		if (arrayTree[i].value != NULL) {        //if value if not null, it's leafNode
			HC[j] = codeArray[i].ptr;        //write into Huffman code
			value[j] = arrayTree[i].value;        //write into leafNode value
			j++;
		}
	}
}


/*for decompressing, reBuild HuffmanTree*/
void HuffmanTree::ReBuildHuffmanTree(const char *str, const int *arr, int i) {
	if (str[strIndex] != '\0') {
		if (str[strIndex] == '1') {                                  //judge if it is childNode
			strIndex++;
			int m = nodeIndex;
			int n = nodeIndex + 1;
			arrayTree[i].lchild = m;                                       //define lchild
			arrayTree[i].rchild = n;                                      //define rchild
			arrayTree[m].parent = arrayTree[n].parent = i;  //define the parent node
			nodeIndex += 2;
			arrayTree[i].value = NULL;                                 //set parentNode value to NULL

			ReBuildHuffmanTree(str, arr, m);
			ReBuildHuffmanTree(str, arr, n);
		}
		else {                                                         //judge it it is leafNode
			strIndex++;
			arrayTree[i].value = arr[leafArrIndex];                     //write leafNode information
			leafArrIndex++;
			return;
		}
	}
	else return;
}



/*decompress file*/
void UncompressFile(const char *geneFilename, const char *codeFilename, const char *backFilename) {
	int Len = 0;                                      //get the number of all nodes
	char treeStruCode[255];  					//read TreeStructure coding from tree information file
	int leafASCII[255] = { 0 };                //read ASCII code of the leaf character from the tree information file
	int i = 0;
	HuffmanTree HT;

	FILE *fr;
	fr = fopen(codeFilename, "rb");
	fscanf(fr, "%d", &Len);                  //read the number of all node into len
	fscanf(fr, "%s", treeStruCode);
	while (fscanf(fr, "%d", &leafASCII[i]) != EOF)    //judge if read all leaf node
		i++;

	int cur = (Len + 1) / 2;
	HT.currsize = cur;
	HT.ReBuildHuffmanTree(treeStruCode, leafASCII, Len - 1);   //root node in the end of arrayTree

	HT.HuffmanCoding();

	char *HC[200] = { 0 };
	char value[200] = { 0 };
	HT.ReBuildArray(HC, value, Len);

	ifstream inFile(geneFilename, ios::in | ios::binary);
	ofstream outFile(backFilename, ios::out | ios::binary);
	if (!inFile || !outFile) {
		cerr << "Open file fail!" << endl;
		exit(1);
	}

	int decToBinInt[10];                 //save,  convert decimal read from compressed file to binary (int type)
	char decToBinChar[1000];       //save,  convert decimal read from compressed file to binary (char type)
	i = 0;
	int j = 0, k, temp = 0;
	while (!inFile.eof()) {
		int data = inFile.get();
		if (data == -1) break;

		memset(decToBinInt, 0, sizeof(decToBinInt));
		i = 0;
		while (data) {                                  //convert compressed decimal to binary
			decToBinInt[i++] = data % 2;
			data = data / 2;
		}

		i = temp;
		for (k = 7; k >= 0; i++, k--) {     //convert binary number into character
			if (decToBinInt[k])
				decToBinChar[i] = '1';
			else
				decToBinChar[i] = '0';

			decToBinChar[i + 1] = '\0';

			for (j = 0; j < HT.currsize; j++) {
				if (strcmp(decToBinChar, HC[j]) == 0) {
					outFile.write((char*)&value[j], sizeof(value[j]));
					i = -1;
				}
			}
		}
		if (i) temp = i;
		else temp = 0;
	}
	inFile.close();
	outFile.close();
	fclose(fr);
	cout << "Decompression Successful!" << endl;
}

int main() {
	int choice;
	char path[100];
	cout << "/********** File compression and decompression using Huffman Code **********/" << endl;
	while (true) {
		cout << "File Compression [1]" << endl;
		cout << "File Decompression [2]" << endl;
		cout << "exit [3]" << endl;
		while (true) {
			cout << "Please enter the number > ";
			cin >> choice;
			if (choice < 1 || choice >3)
				cout << "Invalid Input£°\n" << endl;
			else
				cin.get();
			break;
		}
		switch (choice) {
		case 1:
			cout << "Input the file path to be compressed > ";
			cin.get(path, 100);
								// use your own local path
			CompressFile(path, "/Users/hezhu/Desktop/Huffman/testCompress.bin", "/Users/hezhu/Desktop/Huffman/Compress.txt");
			
			break;
		case 2:
			cout << "Input the file path to be decompressed > ";
			cin.get(path, 100);
			
			UncompressFile(path, "/Users/hezhu/Desktop/Huffman/testCompress.bin", "/Users/hezhu/Desktop/Huffman/Uncompress.txt");
			
			break;
		case 3:
			system("pause");
			return 0;
		default:
			cout << "Wrong file path£°" << endl;
			break;
		}
		cout << endl;
	}
}
