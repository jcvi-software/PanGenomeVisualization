/*              Anika Verma                */
/*   PanGenome Visualization, JCVI 2015    */

#include<iostream>
#include<conio.h>
#include<fstream>
#include<string>
#include<vector>


using namespace std;

int main(){

	string genome;
	string stopcore;
	string ucore;

	cout << "--- PASTE THE CODE FROM trackList.txt IN YOUR tracklist.json FILE ---" << endl << endl;
	cout << "Which Genome file is this? (ex. E00001)" << endl;
	getline(cin, genome);
	cout << "What color would you like the STOP_CORE subfeature to be? (ex. dblhelix, est, transcript-cds )" << endl;
	getline(cin, stopcore);
	cout << "What color would you like the U_CORE subfeature to be? (ex. dblhelix, est, transcript-cds ) " << endl;
	getline(cin, ucore);

	ofstream file("genomeTrackList.txt");

	file << "{" << endl;
	file << " \"style\" : {  " << endl;
	file << " \t \"className\" : \"STOP_CORE, U_CORE, Break\", " << endl;
	file << " \t \"subfeatureClass\" : { " << endl;
	file << " \t \t \"STOP_CORE\" : \""<< stopcore <<"\", " << endl;
	file << " \t \t \"U_CORE\" : \"" << ucore << "\", " << endl;
	file << "\t }" << endl;
	file << " }, " << endl;
	file << " \"key\" : \"" << genome << " Genome Track\", " << endl;
	file << " \"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", " << endl;
	file << " \"trackType\" : \"null\", " << endl;
	file << " \"urlTemplate\" : \"tracks/" << genome << "Track/{refseq}/trackData.json\", " << endl;
	file << " \"compress\" : \"0\", " << endl;
	file << " \"label\" : \"" << genome << " Track\", " << endl;
	file << " \"type\" : \"CanvasFeatures\", " << endl;
	file << " \"category\" : \"" << genome << "\" " << endl;
	file << "}," << endl;

	cout << "Complete!  Check your genomeTrackList.txt file for the results. " << endl;
	_getch();
	return 0;
}
