#include"HuffmanCompress_Un.h"

void HuffmanTree::Process(const char *inFilename, const char *outFilename)
{
	ifstream inFile(inFilename, ios::in | ios::binary);       
	if (!inFile) {
		cerr << " Can not open!" << endl;
		exit(1);
	}
	char ch = inFile.get();   

	while (!inFile.eof()) {           
		int pos = FindPosition(ch);
		if (pos != -1)
			arrayTree[pos].weight++;
		else
			Insert(ch, 1);        

		ch = inFile.get();
	}
	BuildHuffmanTree();                   
	HuffmanCoding();                 
	WwriteCodeToFile(outFilename);
	inFile.close();
}



int HuffmanTree::FindPosition(const char &ch) const {
	for (int i = 0; i < currsize; i++)
		if (arrayTree[i].value == ch)
			return i;
	return -1;
}



void HuffmanTree::Insert(const char &data, const int &wt) {
	if (2 * currsize - 1 > maxSize) return;
	arrayTree[currsize].value = data;
	arrayTree[currsize].weight = wt;
	currsize++;
}



void HuffmanTree::BuildHuffmanTree() {
	int i = currsize;
	int wt1, wt2;
	int s1 = 0, s2 = 0;
	while (i < 2 * currsize - 1) {
		wt1 = wt2 = numeric_limits<int>::max();     
		for (int k = 0; k < i; k++) {
			if (arrayTree[k].weight < wt1 && arrayTree[k].parent == -1) {
				wt2 = wt1;
				s2 = s1;
				wt1 = arrayTree[k].weight;                      
				s1 = k;
			}
			else if (arrayTree[k].weight < wt2 && arrayTree[k].parent == -1) {
				wt2 = arrayTree[k].weight;                      
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


void HuffmanTree::HuffmanCoding()
{
	char tempCodeArray[100] = { 0 };                       
	if (codeArray)
	{
		delete[]codeArray;
		codeArray = NULL;
	}
	codeArray = new HuffmanCode[2 * currsize - 1];  
	if (!codeArray)
		return;

	int Len = 2 * currsize - 1;                        
	int p, tmp_index;
	for (int i = 0; i < Len; i++) {
		tmp_index = i;
		p = arrayTree[i].parent;
		int n = 0;
		while (p != -1 && p < Len) {
			if (arrayTree[p].lchild == tmp_index)
				tempCodeArray[n++] = '0';        
			else if (arrayTree[p].rchild == tmp_index)
				tempCodeArray[n++] = '1';         

			tmp_index = p;
			p = arrayTree[tmp_index].parent;
		}
		for (int t = 0; t < n; t++) {                        
			codeArray[i].ptr[t] = tempCodeArray[n - 1 - t];

		}
		codeArray[i].ptr[n] = '\0';
	}
}



void HuffmanTree::WwriteCodeToFile(const char *outFilename) {


	ofstream fw(outFilename, ios::out | ios::binary);
	if (fw)
	{
		fw << currsize * 2 - 1;                  
		WriteTreeStructureToFile(currsize * 2 - 2, fw);
		fw << " ";
		WriteLeafToFile(currsize * 2 - 2, fw);               
		fw.close();
	}

}



void HuffmanTree::WriteTreeStructureToFile(const int &num, ofstream &fw) {
	if (arrayTree[num].lchild == -1) {
		fw << "0";                        
		return;
	}
	else
		fw << "1";                        

	WriteTreeStructureToFile(arrayTree[num].lchild, fw);
	WriteTreeStructureToFile(arrayTree[num].rchild, fw);
}




void HuffmanTree::WriteLeafToFile(const int &num, ofstream&fw)
{
	if (arrayTree[num].lchild == -1) {          
		fw << arrayTree[num].value << " "; 
		return;
	}
	WriteLeafToFile(arrayTree[num].lchild, fw);
	WriteLeafToFile(arrayTree[num].rchild, fw);
}


void CompressFile(const char *sourceFilename, const char *codeFilename, const char *geneFilename) {
	HuffmanTree tree;
	tree.Process(sourceFilename, codeFilename);

	ifstream inFile(sourceFilename, ios::in | ios::binary);
	ofstream outFile(geneFilename, ios::out | ios::binary);
	if (!inFile || !outFile) {
		cerr << "Open file fail!" << endl;
		exit(1);
	}

	int decimal = 0;                        
	int i, j, k = 0;
	int bits = 0;                                
	int temp[10000];
	memset(temp, 0, sizeof(temp));


	while (!inFile.eof()) {
		decimal = 0;
		unsigned char one_byte = inFile.get();
		for (i = 0; i < tree.currsize; i++) {
			if (one_byte == tree.arrayTree[i].value) {
				bits += strlen(tree.codeArray[i].ptr);
				int len = strlen(tree.codeArray[i].ptr);  
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
						temp[j] = 0;                  

					for (j = 8; j < k; j++)
						temp[j - 8] = temp[j];        
					k = bits = j - 8;

					unsigned char c = decimal;   
					outFile.write((char*)&c, sizeof(c));
				}
			}
		}
	}
	if (bits) {             
		decimal += temp[0] * 128 + temp[1] * 64 + temp[2] * 32 + temp[3] * 16 + temp[4] * 8 + temp[5] * 4 + temp[6] * 2 + temp[7] * 1;
		unsigned char c = decimal;
		outFile.write((char*)&c, sizeof(c));
	}
	inFile.close();
	outFile.close();
	cout << "Compression Successful!" << endl;
}




void HuffmanTree::ReBuildArray(char **HC, char *value, int Num) {
	int j = 0;
	for (int i = 0; i < Num; i++) {
		if (arrayTree[i].value != NULL) {        
			HC[j] = codeArray[i].ptr;        
			value[j] = arrayTree[i].value;        
			j++;
		}
	}
}



void HuffmanTree::ReBuildHuffmanTree(const char *str, const int *arr, int i) {
	if (str[strIndex] != '\0') {
		if (str[strIndex] == '1') {                                           
			strIndex++;
			int m = nodeIndex;
			int n = nodeIndex + 1;
			arrayTree[i].lchild = m;                                       
			arrayTree[i].rchild = n;                                        
			arrayTree[m].parent = arrayTree[n].parent = i;  
			nodeIndex += 2;
			arrayTree[i].value = NULL;                                        
			ReBuildHuffmanTree(str, arr, m);
			ReBuildHuffmanTree(str, arr, n);
		}
		else {                                                                     
			strIndex++;
			arrayTree[i].value = arr[leafArrIndex];                    
			leafArrIndex++;
			return;
		}
	}
	else return;
}




void UncompressFile(const char *geneFilename, const char *codeFilename, const char *backFilename) {
	int Len = 0;                                      
	char treeStruCode[255];  
	int leafASCII[255] = { 0 };                
	int i = 0;
	HuffmanTree HT;

	FILE *fr;
	fr = fopen(codeFilename, "rb");
	fscanf(fr, "%d", &Len);                  
	fscanf(fr, "%s", treeStruCode);
	while (fscanf(fr, "%d", &leafASCII[i]) != EOF)    
		i++;

	int cur = (Len + 1) / 2;
	HT.currsize = cur;
	HT.ReBuildHuffmanTree(treeStruCode, leafASCII, Len - 1);   

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

	int decToBinInt[10];              
	char decToBinChar[1000];       
	i = 0;
	int j = 0, k, temp = 0;
	while (!inFile.eof()) {
		int data = inFile.get();
		if (data == -1) break;

		memset(decToBinInt, 0, sizeof(decToBinInt));
		i = 0;
		while (data) {                                  
			decToBinInt[i++] = data % 2;
			data = data / 2;
		}

		i = temp;
		for (k = 7; k >= 0; i++, k--) {     
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
				cout << "Invalid Input!\n" << endl;
			else
				cin.get();
			break;
		}
		switch (choice) {
		case 1:
			cout << "Input the file path to be compressed > ";
			cin.get(path, 100);
			CompressFile(path, "/Users/hezhu/Desktop/For_Job/Pre_Autumn2020/changshi/testCompress.bin", "/Users/hezhu/Desktop/For_Job/Pre_Autumn2020/changshi/Compress.txt");
			break;
		case 2:
			cout << "Input the file path to be decompressed > ";
			cin.get(path, 100);
			UncompressFile(path, "/Users/hezhu/Desktop/For_Job/Pre_Autumn2020/changshi/testCompress.bin", "/Users/hezhu/Desktop/For_Job/Pre_Autumn2020/changshi/Uncompress.txt");
			break;
		case 3:
			system("pause");
			return 0;
		default:
			cout << "Wrong file path!" << endl;
			break;
		}
		cout << endl;
	}
}