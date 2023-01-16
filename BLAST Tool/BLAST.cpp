/*
Author: Rajiv Snape

This is project is a BLAST tool. BLAST stands for Basic Local Alignment Search Tool is one of the most important tools in bioinformatics.
The tool works by entering a DNA sequence(s) (or protein if option is avialabe) and then search a database of thousands of nucleotides. It will then
report if sequence(s) are found or not. If so then their locations in the database will be reported. Before this is done the user provided sequences and te enter
bases are converted to K-mers of user determined length.

*/

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <map>

using namespace std;

void read_in_database_file(string db_file);
void finding_db_kmers(vector<string> database);
int kmer_length();
string read_in_query();
void finding_query_kmers(multimap<string, pair<size_t, size_t>> kmer_db_map, size_t kmer_size);
void kmer_query_lookup(multimap<string, pair<size_t, size_t>> kmer_db_map, size_t kmer_size, vector <string> query_kmers);

vector<string> database;
string query_seq = "";


int kmer_input() 
{// asks the user for the lenght of k-mer 
	int kmer_length;
	cout << "Please enter k-mer length" << endl;
	cin >> kmer_length;

	cout << endl;

	return kmer_length;
}

string read_in_query() 
{
	// asks the user for the name of query file containing the desired seqeuence they want to locate in the database
	string query_file = "";
	cout << "Please enter the name of your query file you want to search for: " << endl;

	getline(cin >> ws, query_file);

	ifstream fin;
        fin.open(query_file);
        if (!fin.is_open())
        {//if
                cerr << "Error did not open file" << endl;
                exit(1);
        }//if

	string line = "";
	string header = "";
	string sequence = "";


    while(getline(fin,line))
	{
        if( line.empty() || line[0] == '>' )
        { // Identifier marker
            if(!header.empty() )
            { // Print out what we read from the last entry
                header.clear();
                query_seq += sequence;
            }

            if( !line.empty() )
            {
                header = line.substr(1);
            }

            sequence.clear();
        }

        else if(!header.empty())
        {
            line = line.substr(0, line.length() -1);
            if(line.find(' ') != string::npos )
            { // Invalid sequence--no spaces allowed
                header.clear();
                sequence.clear();
            }

            else
            {
                sequence += line;
            }
        }
    }

    if(!header.empty() )
    { // Print out what we read from the last entry
        query_seq += sequence;
    }

	return query_seq;
}


void kmer_query_lookup(multimap<string, pair<size_t, size_t>> kmer_db_map, size_t kmer_size, vector <string> query_kmers) {

	map<size_t, pair<size_t, size_t>> database_num_counter_index; // hash map that keeps track of databases

    for (size_t d_num = 0; d_num < 13; d_num++)
    { // loops through index of the database
		database_num_counter_index.insert({d_num, make_pair(0,0)});
    }

	for (auto i : query_kmers)
	{ // Loops throughs all the kmers creadted from the query
		auto it = kmer_db_map.equal_range(i); // using the properties of a hash/multi-map the query k-mers are found and located in the multi-map of k-mers

		for (auto itr = it.first; itr != it.second; itr++)
		{
			pair<size_t, size_t> start_and_counter = database_num_counter_index[itr->second.first];

			if (start_and_counter.first == 0)
			{
				start_and_counter.first = itr->second.second;
			}
			else if ((start_and_counter.first + start_and_counter.second) == itr->second.second - 1)
			{
				start_and_counter.second =  start_and_counter.second + 1;
			}
			else if (start_and_counter.second > 0)
			{
				string super_kmer;

				for (size_t x = start_and_counter.first; x < start_and_counter.first +  start_and_counter.second + kmer_size; x++)
				{
					super_kmer += database[itr->second.first][x];
				}

				cout << "The following sequence of your query : " <<
				"\n" << super_kmer << "\n" <<
				"has been detected at Database " << itr->second.first+1 <<
				" " << start_and_counter.first+1 << " - " << start_and_counter.first +  start_and_counter.second + kmer_size <<
				" " << "|| Query " << query_seq.find(super_kmer)+1 << " - " << query_seq.find(super_kmer)+super_kmer.length() << endl;

				string concensus_db = database[itr->second.first].substr(start_and_counter.first, super_kmer.length());
				string concensus_q = query_seq.substr(query_seq.find(super_kmer), super_kmer.length());
				string concensus = "";

				for (size_t db = 0; db <= concensus_db.length(); db++)
				{
					char index1 = concensus_db[db];
					char index2 = concensus_q[db];

					if (index1 == index2)
					{
						concensus += index1;
					}
					else
					{
						concensus = '.';
					}
				}

				cout << "Matching database sequence " << "\n" << concensus_db << endl;
				cout << "Concensus: " << "\n" << concensus << endl;

				cout << endl;
                cout << "=================================================================================" << endl;
                cout << endl;

				start_and_counter.first = itr->second.second;
				start_and_counter.second = 0;
			}
				database_num_counter_index[itr->second.first] = start_and_counter;
		}
	}

	for (size_t d_num = 0; d_num < 13; d_num++)
    {
		pair<size_t, size_t> start_and_counter = database_num_counter_index[d_num];

		if (start_and_counter.second > 0)
        {

			string super_kmer;

			for (size_t x = start_and_counter.first; x < start_and_counter.first +  start_and_counter.second + kmer_size; x++)
            {
                super_kmer += database[d_num][x];
            }

            cout << "The following sequence of your query : " <<
			"\n" << super_kmer << "\n" << "has been detected at database " << d_num+1 <<
            " " << start_and_counter.first+1 << " - " << start_and_counter.first +  start_and_counter.second + kmer_size <<
			" " << "|| Query " << query_seq.find(super_kmer) << " - " << query_seq.find(super_kmer)+super_kmer.length() << endl;

			string concensus_db = database[d_num].substr(start_and_counter.first, super_kmer.length());
            string concensus_q = query_seq.substr(query_seq.find(super_kmer), super_kmer.length());
            string concensus = "";

            cout << endl;

            for (size_t db = 0; db <= concensus_db.length(); db++)
            {
                char index1 = concensus_db[db];
                char index2 = concensus_q[db];

                if (index1 == index2)
                {
                    concensus += index1;
                }
                else
                {
                    concensus = '.';
                }
            }

            cout << "Matching database sequence: " << "\n" << concensus_db << endl;
            cout << "Concensus: " << "\n" << concensus << endl;

            cout << endl;
			cout << "=================================================================================" << endl;
			cout << endl;
		}

	}
}

void finding_query_kmers(multimap<string, pair<size_t, size_t>> kmer_db_map, size_t kmer_size) 
{

	string query_seq = read_in_query(); // read in and storing user entered query sequences that will be searched for in the database
	cout << "Your query sequence is " << query_seq << endl;
	cout << endl;

	vector <string> query_kmers;
	for (size_t i = 0; i <= query_seq.length(); i++) // looping through each character in the query sequence
	{
		string temp = query_seq.substr(i, kmer_size); // sliding window of kmer soze
		if (temp.length() == kmer_size) // makes sure the window never runs less than the given kmer size
		{
			query_kmers.push_back(temp); // adds generated kmers to vector where they are stored.
		}
	}

	kmer_query_lookup(kmer_db_map, kmer_size, query_kmers);
}

void finding_db_kmers(vector<string> database) 
{ // This function turns the entire database into K-mers of the user size

	multimap<string, pair<size_t, size_t>> kmer_db_map; // hash table for db kmers. Key = kmer, value database location & position in database

	size_t kmer_size = kmer_input(); // receives k-mer siae from user

	for (size_t y = 0; y < database.size(); y++)
	{// loops through the entire database
		for (size_t i = 0; i <= database[y].length(); i++) // loops though each database sequence
		{
			string kseq = database[y].substr(i, kmer_size); // sliding window sequence
			pair<size_t, size_t> location; 
			location = make_pair(y,i); // create a pair storing the k-mer in the database and sequene

			if (kseq.length() == kmer_size)
			{// to ensure the window length didn't isn't less than the k-mer size
				kmer_db_map.insert(pair<string, pair<size_t, size_t>>(kseq, location));

			}

		}
	}

	finding_query_kmers(kmer_db_map, kmer_size);
}

void read_in_database_file(string db_file) {

	ifstream fin;
        fin.open(db_file);
        if (!fin.is_open())
        {//if
                cerr << "Error did not open file" << endl;
                exit(1);
        }//if

	string line = "";

	string seqnames = "";
	string sequence = "";

	while(getline(fin,line))
	{
       	if(line.empty() || line[0] == '>')
		{ // Identifier marker
            if(!seqnames.empty() )
			{ // Print out what we read from the last entry
                seqnames.clear();
				database.push_back(sequence);
            }
 
			if( !line.empty() )
			{
            	seqnames = line.substr(1);
            }
            sequence.clear();
		}

		else if(!seqnames.empty()) 
		{
			line = line.substr(0, line.length() -1);
           	if(line.find(' ') != string::npos )
			{ // Invalid sequence--no spaces allowed
				seqnames.clear();
                sequence.clear();
        	}

			else
			{
			    sequence += line;
            } 
		}
    }

    if(!seqnames.empty())
	{ // Print out what we read from the last entry
		database.push_back(sequence);
	}

	finding_db_kmers(database);
}

int main() {

	cout << "\n" << endl;

	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "------------------------------------ BLAST -----------------------------------" << endl;
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

	cout << endl;

	cout << "Welcome to BLAST" << endl;
	cout << "Please enter the name of your database file" << endl;
	string db_file = "";
	getline(cin, db_file);

	cout << endl;

	read_in_database_file(db_file);
}
