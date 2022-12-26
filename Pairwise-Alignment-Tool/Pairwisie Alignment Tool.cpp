/*
Author: Rajiv Snape

This is project is a pairwise alignment tool for DNA sequences and Protein sequences.
At the start of the program you choose wether you want to do a DNA alignment or protein alignment.
For either case the user will enter a FASTA file (or a file in FASTA format) with only 2 sequences.
Then algoirthm will then perform a Needleman – Wunsch global alignment. For DNA alignment, the user enters mismatch
gap, and match scores. For protein alignment user enters the gap score while the othe values will be derived from
a user entered PAM scroing matrix.

*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <climits>

using namespace std;

void input_DNA();
void input_Protein();
void read_in_DNA(string);
void read_in_Protein(string);
void Pairwise_aligning_nuc(string Seq1, string Seq2, int Seq1_length, int Seq2_length);
void Pairwise_alignment_prot(string Prot1_Seq, string Prot2_Seq, int Prot1_len, int Prot2_len);

void Pairwise_aligning_nuc(string Seq1, string Seq2, int Seq1_length, int Seq2_length){ // doing the DNA alignment

	int gap_pen{};
	int match_score{};
	int mismatch{};


	// Asking the user for mismatchm gap, and match scores
	cout << "\nPlease enter gap penalty" << endl;
   	cin >> gap_pen;
   	cout << "Please enter score for a match" << endl;
   	cin >> match_score;
   	cout << "Please enter penalty for mismatch" << endl;
   	cin >> mismatch;


	// Declartion of scoring matrix in which alignment will be performed.
   	int matrix[Seq1_length+1][Seq2_length+1] {0};

	// These two for loops apply the gap penalty along the borders of the matrix as per Needleman–Wunsch algorithm
   	for (int i {1}; i <= Seq1_length; i++)
   	{
		matrix[i][0] = gap_pen  * i;
   	}
    	for (int j {1}; j <= Seq2_length; j++)
	{
       		matrix[0][j] = gap_pen * j;
   	}

	// These two for loops perform the calculatiions. They loop through each provided sequence at each base is compared.  
	for(int i {0}; i <= Seq1_length; i++)
    	{
        	for(int j {0}; j <= Seq2_length; j++)
        	{
            		if (i > 0 && j > 0)
            		{
				// If the sequences share the same base at certain location the corresponding coordinate in the matrix is filled with match score + the diagonal value.
                		if (Seq1[i-1] == Seq2[j-1])
    		    		{
                   			matrix[i][j] = matrix[i-1][j-1] + match_score;
                		}
				// If not a match the coordinate of the matrix will be filled with highest value of either horizontal gap score, vertical gap score or
				// diagonal value + mismatch
               			else
                		{
					// Highest value is chosen out of mismtach or gap penalties
                    			matrix[i][j] = max({matrix[i-1][j-1] + mismatch  ,
                                                matrix[i-1][j] + gap_pen ,
                                                matrix[i][j-1] + gap_pen});
                		}
            		}
        	}
    	}
	cout <<  endl;

	// This finds the score of the alignment which will be the value at the most bottom right
	int terminal = matrix[Seq1_length+1-1][Seq2_length+1-1];
	cout << terminal << endl;

	int temp = terminal;
	int decrement_X = Seq1_length+1-1;
	int decrement_Y = Seq2_length+1-1;

	string Seq1_aligned = "";
	string Seq2_aligned = "";

	int inc_seq1 = Seq1_length-1;
	int inc_seq2 = Seq2_length-1;

	// starting from the bottom right cell we back calculate where the current cell comes from
	while (!(decrement_X == 0 || decrement_Y == 0))
	{
		// If the current cell's value is equal to the cell above it + the gap penalty, a gap is introduced to the top sequence, the bottom remains the same
		if (terminal == matrix[decrement_X][decrement_Y-1] + gap_pen) 
		{
            		temp = matrix[decrement_X][decrement_Y-1];
            		terminal = temp;
            		decrement_Y--;
			Seq1_aligned += '_';
			Seq2_aligned += Seq2[inc_seq2];
			inc_seq2--;
        	}

		// If the current cell's value is equal to the cell to the left + the gap penalty, a gap is introduced to the bottom sequence, the top remains the same
		else if (terminal == matrix[decrement_X-1][decrement_Y] + gap_pen) 
		{
            		temp = matrix[decrement_X-1][decrement_Y];
            		terminal = temp;
            		decrement_X--;
            		Seq1_aligned += Seq1[inc_seq1];
            		Seq2_aligned += '_';
			inc_seq1--;
        	}

		// If the current cell's value is equal to the diagonal cell's value + the match score, The nucleotide bases of both sequences are added
		else if (terminal == matrix[decrement_X-1][decrement_Y-1] + match_score) 
		{
			temp = matrix[decrement_X-1][decrement_Y-1];
			terminal = temp;
			decrement_X--;
			decrement_Y--;
			Seq1_aligned += Seq1[inc_seq1];
			Seq2_aligned += Seq2[inc_seq2];
			inc_seq1--;
			inc_seq2--;
		}

		// If the current cell's value is equal to the diagonal cell's value + the mismatch score, a gap is introduced to both sequences
		else if (terminal == matrix[decrement_X-1][decrement_Y-1] + mismatch) 
		{
            		temp = matrix[decrement_X-1][decrement_Y-1];
            		terminal = temp;
			decrement_X--;
            		decrement_Y--;
			Seq1_aligned += Seq1[inc_seq1];
            		Seq2_aligned += Seq2[inc_seq2];
			inc_seq1--;
            		inc_seq2--;
        	}

	}

	// These two while loops handle to condition of gaps being at the start of the aligned sequences
	// If the current matrix corrdinate has a y value of zero introduce gaps to the bottom sequence and added the base of the current top sequence
	while (decrement_X >= 1 && decrement_Y < 1)
    	{
    		Seq1_aligned += Seq1[inc_seq1];
		Seq2_aligned += '_';
		decrement_X--;
    	}

	// If the current matrix corrdinate has a x value of zero introduce gaps to the top sequence and added the base of the current bottom sequence
   	 while (decrement_Y >= 1 && decrement_X < 1)
    	{
		Seq1_aligned += '_';
		Seq2_aligned += Seq2[inc_seq2];;
		decrement_Y--;
    	}

	// These reverse the order of aligned sequences.
	reverse(Seq1_aligned.begin(), Seq1_aligned.end());
	reverse(Seq2_aligned.begin(), Seq2_aligned.end());

	cout << "Performing alignment" << " .................................................. 100%" << endl;
	cout << "Alignment complete! Your aligned sequences are:" << endl;

	cout << "Sequence 1: " << Seq1_aligned << endl;
    	cout << "Sequence 2: " << Seq2_aligned << endl;
} // doing the alignment


void Pairwise_alignment_prot(string Prot1_Seq, string Prot2_Seq, int Prot1_len, int Prot2_len) 
{// Performing protein alignment

	int gap_pen {};

	// Asking the user for gap penalty which is then stored in the gap_pen variable
	cout << "\nPlease enter a gap penalty: " << endl;
	cin >> gap_pen;
	cout << endl;

	// Creating and initializing the alignment matrix based on the lengths of provided protein sequences
	int matrix[Prot1_len+1][Prot2_len+1] {0};

	// Applying the gap penalty along the edges as per the Needleman–Wunsch algorithm with these two for loops
    for (int i {1}; i <= Prot1_len; i++)
    {
        matrix[i][0] = gap_pen  * i;
    }
    for (int j {1}; j <= Prot2_len; j++)
    {
        matrix[0][j] = gap_pen * j;
    }

	// The following statements reads in a PAM matrix file that is provdie by the user which conatins the scores for all the possible combination between amino acids
	// File is requested and stored in PAM_file variable
	string PAM_file;
	cout << "Enter the name of your PAM matrix file " << endl;
	cin >> PAM_file;
	fstream myfile;
    myfile.open(PAM_file, ios::in);

	if (!myfile.is_open())
    {
        cerr << "Error can not open file" << endl;
        exit(1);
    }

	// pam_matrix is declared
	const int len = 20;
    int pam_matrix[len][len];

    string line;

	// getting rid of the first line of pam matrix file, by using getline to grab it
	// The first line constains amino acid labels. The goal is to store only the values of table into matrix created above
    getline(myfile,line); 

    for (int i = 0; i < len; i++)
    {
        char Row_first_letter;
        myfile >> Row_first_letter; // using the >> to grab the first char of each row leaving the pam matrix with just numbers

        for(int j = 0; j < len; j++)
        {
            myfile >> pam_matrix[i][j]; // store file contents into pam_matrix variable
        }
    }

	cout << endl;

	int match_score {};
	int mismatch {};


	// Amino acid labels are stored in 2 arrays the same order they would appear in the table
	char row_labels [20] {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
	char col_labels [20] {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};


	// These two for loops perform the calculatiions of the alignment matrix. They loop through each provided sequence at each base is compared.
	for(int i {0}; i <= Prot1_len; i++)
    {
        for(int j {0}; j <= Prot2_len; j++)
        {
            if (i > 0 && j > 0)
            {
				// If two amino acids at a particular position is identiical their match value is added to the of the diagonal cell
                if (Prot1_Seq[i-1] == Prot2_Seq[j-1])
                {
					// the row and col variables were created to store the index of an amino acid in the labels array. This index would correspond to the index of the
					// amino acid in the header row/colunm of the PAM matrix file
					int row = distance(row_labels, find(row_labels, row_labels + 20, Prot1_Seq[i-1]));
                    int col = distance(col_labels, find(col_labels, col_labels + 20, Prot2_Seq[j-1]));

					// the match value is then the combined index of both amino acid labels in the pam_matrix 2D-array data structure craeted
					match_score = pam_matrix[row][col];

					// If the sequences share the same base at certain location the corresponding coordinate in the matrix is filled with match score + the 
					// diagonal value.
					matrix[i][j] = matrix[i-1][j-1] + match_score;
                }

				// If not a match the coordinate of the matrix will be filled with highest value of either horizontal gap score, vertical gap score or
				// diagonal value + mismatch
				else if (Prot1_Seq[i-1] != Prot2_Seq[j-1])
				{

					// the row and col variables were created to store the index of an amino acid in the labels array. This index would correspond to the index of the
					// amino acid in the header row/colunm of the PAM matrix file
					int row = distance(row_labels, find(row_labels, row_labels + 20, Prot1_Seq[i-1]));
					int col = distance(col_labels, find(col_labels, col_labels + 20, Prot2_Seq[j-1]));

					// the mismatch value is then the combined index of both amino acid labels in the pam_matrix 2D-array data structure craeted
					mismatch = pam_matrix[row][col];

					// Highest value is chosen out of mismtach or gap penalties
					matrix[i][j] = max({matrix[i-1][j-1] + mismatch,
                                                matrix[i-1][j] + gap_pen ,
                                                matrix[i][j-1] + gap_pen});
				}
            }
        }
    }

		// This finds the score of the alignment which will be the value at the most bottom right
        int terminal = matrix[Prot1_len+1-1][Prot2_len+1-1]; // rows and colums are length+1 long. and to get the edge value it's col/row length -1.

        int temp = terminal;
        int decrement_X = Prot1_len+1-1;
        int decrement_Y = Prot2_len+1-1;

        string Prot1_aligned = "";
        string Prot2_aligned = "";

        int inc_seq1 = Prot1_len-1;
        int inc_seq2 = Prot2_len-1;

		// starting from the bottom right cell we back calculate where the current cell comes from
        while (!(decrement_X == 0 || decrement_Y == 0))
        {//while

			int match_index_1;
        	int match_index_2;

			// Since the match and mismatch values come from the PAM matrix, these values need to be changed and update for every base pair combination
			// As the program backtracks it changes and retreives the combination score of each amino acid pair in the PAM matrix
			// which is sotred as match ot mismatch
			for (int i = 0; i < 20; i++)
            {
               	if (row_labels[i] == Prot1_Seq[inc_seq1])
                {
                    match_index_1 = i;
                }
            }
            for (int i = 0; i < 20; i++)
            {
                if (col_labels[i] == Prot2_Seq[inc_seq2])
                {
                    match_index_2 = i;
                }
            }
            match_score = pam_matrix[match_index_1][match_index_2];

			int mismatch_index_1;
			int mismatch_index_2;

			for (int i = 0; i < 20; i++)
			{
				if(row_labels[i] == Prot1_Seq[inc_seq1])
				{
					mismatch_index_1 = i;
				}
			}
			for (int i = 0; i < 20; i++)
			{
				if(col_labels[i] == Prot2_Seq[inc_seq2])
				{
					mismatch_index_2 = i;
				}
			}
			mismatch = pam_matrix[mismatch_index_1][mismatch_index_2];

			// If the current cell's value is equal to the cell above it + the gap penalty, a gap is introduced to the top sequence, the bottom remains the same
            if (terminal == matrix[decrement_X][decrement_Y-1] + gap_pen) 
			{
                temp = matrix[decrement_X][decrement_Y-1];
                terminal = temp;
                decrement_Y--;
                Prot1_aligned += "_";
                Prot2_aligned += Prot2_Seq[inc_seq2];
                inc_seq2--;
            }

			// If the current cell's value is equal to the cell to the left + the gap penalty, a gap is introduced to the bottom sequence, the top remains the same
            else if (terminal == matrix[decrement_X-1][decrement_Y] + gap_pen) 
			{
                temp = matrix[decrement_X-1][decrement_Y];
                terminal = temp;
                decrement_X--;
                Prot1_aligned += Prot1_Seq[inc_seq1];
                Prot2_aligned += "_";
                inc_seq1--;
            }

			// If the current cell's value is equal to the diagonal cell's value + the match score, The nucleotide bases of both sequences are added
            else if (terminal == matrix[decrement_X-1][decrement_Y-1]+ match_score) 
			{
                temp = matrix[decrement_X-1][decrement_Y-1];
                terminal = temp;
                decrement_X--;
                decrement_Y--;
            	Prot1_aligned += Prot1_Seq[inc_seq1];
                Prot2_aligned += Prot2_Seq[inc_seq2];
				inc_seq1--;
                inc_seq2--;
            }

			// If the current cell's value is equal to the diagonal cell's value + the mismatch score, a gap is introduced to both sequences
        	else if (terminal == matrix[decrement_X-1][decrement_Y-1] + mismatch) 
			{
                temp = matrix[decrement_X-1][decrement_Y-1];
                terminal = temp;
                decrement_X--;
                decrement_Y--;
                Prot1_aligned += Prot1_Seq[inc_seq1];
                Prot2_aligned += Prot2_Seq[inc_seq2];
                inc_seq1--;
                inc_seq2--;
            }

        } //while

	// These reverse the order of aligned sequences.
	reverse(Prot1_aligned.begin(), Prot1_aligned.end());
    reverse(Prot2_aligned.begin(), Prot2_aligned.end());

    cout << "Performing alignment" << " .................................................. 100%" << endl;
    cout << "Alignment complete! Your aligned sequences are:" << endl;

	cout << endl;
    cout << "Sequence 1: " << Prot1_aligned << endl;
    cout << "Sequence 2: " << Prot2_aligned << endl;

}

void read_in_DNA(string DNA_filename) 
{ // read in the sequences
    fstream myfile;
    myfile.open(DNA_filename, ios::in);

    if (!myfile.is_open())
    {
        cerr << "Error can not open file" << endl;
        exit(1);
    }

	string Seq1{};
    string Seq2{};
    string Seqnames{};

    string temp{};

	string line{};
    while(getline(myfile,line))
	{
        if (line[0] == '>')
		{//if
			Seqnames.append(line);
			if (temp != "")
			{//if
				temp = temp.substr(0,temp.length()-1);
				Seq1.append(temp);
				temp = "";
			}//if
		}//if
		else
		{//else
			temp += line;
		}//else
    	}
    	Seq2.append(temp);

    	int Seq1_length{};
    	int Seq2_length{};

        Seq1_length =  Seq1.length();

        Seq2_length = Seq2.length();

    Pairwise_aligning_nuc(Seq1, Seq2, Seq1_length, Seq2_length);
}


void read_in_Protein(string Protein_filename) 
{ // read in the sequences
    fstream myfile;
    myfile.open(Protein_filename, ios::in);

    if (!myfile.is_open())
    {
        cerr << "Error can not open file" << endl;
        exit(1);
    }

    string Protein_Sequences{};
    string Protein_Seq_names{};

	string Prot1_Seq{};
	string Prot2_Seq{};

    string line{};

	while(getline(myfile,line))
	{
    	if( line.empty() || line[0] == '>' )
		{ // Identifier marker
            if(!Protein_Seq_names.empty() )
			{ // Print out what we read from the last entry
                Protein_Seq_names.clear();
				Prot1_Seq = Protein_Sequences;
            }

            if( !line.empty() )
			{
				Protein_Seq_names = line.substr(1);
            }
       		Protein_Sequences.clear();
		}
		else if(!Protein_Seq_names.empty()) 
		{
			line = line.substr(0, line.length() -1);
           	if(line.find(' ') != string::npos )
			{ // Invalid sequence--no spaces allowed
        		Protein_Seq_names.clear();
				Protein_Sequences.clear();
            }
			else
			{
				Protein_Sequences += line;
            }
        }	
	}

    if(!Protein_Seq_names.empty() )
	{ // Print out what we read from the last entry
		Prot2_Seq = Protein_Sequences;
	}

	cout << "\nSequence 1: " << Prot1_Seq << endl;
	int Prot1_len = Prot1_Seq.length();

	cout << "\nSequence 2: " << Prot2_Seq << endl;
	int Prot2_len = Prot2_Seq.length();

	Pairwise_alignment_prot(Prot1_Seq, Prot2_Seq, Prot1_len, Prot2_len);
}

void input_DNA() {

	string DNA_filename{};

	cout << "\n++++++++++++++++++++++++++++ Pairwise Nucletide Alignment +++++++++++++++++++++++++++++++++" << endl;
    	cout << "\nPlease enter the name of sequence file: " << endl;
    	cin >> DNA_filename;
    	cout << endl;
    	read_in_DNA(DNA_filename);
}

void input_Protein() {

	string Protein_filename{};

        cout << "\n++++++++++++++++++++++++++++ Pairwise Protein Alignment +++++++++++++++++++++++++++++++++" << endl;
        cout << "\nPlease enter the name of sequence file: " << endl;
	cin >> Protein_filename;
	cout << endl;
	read_in_Protein(Protein_filename);
}

int main() {

	char align_choice {};

	do
	{
		cout << "\n=============================================================================" << endl;
		cout << "************************ Raj's Bioinformatics Alignment Tool *****************" << endl;
		cout << "==============================================================================" << endl;

		cout << "\nWould you like to perform a DNA(Nucleotide) Alignment or a Protein(Amino Acid) Alignment?" << endl;
		cout << "Choose N for DNA or P for protein" << endl;
		cin >> align_choice;

		if (align_choice == 'N' || align_choice == 'n')
		{
			input_DNA();
		}

		else if (align_choice == 'P' || align_choice == 'p')
		{
			input_Protein();
		}

		else
		{
			cout << align_choice << " is not a valid choice..... Try again" << endl;
		}

	} while((align_choice != 'N' && align_choice != 'n') && (align_choice != 'P' && align_choice != 'p'));

    	cout << endl;
    	return 0;
}
