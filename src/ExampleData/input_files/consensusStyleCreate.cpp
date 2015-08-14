/*              Anika Verma                */
/*   PanGenome Visualization, JCVI 2015    */

#include<iostream>
#include<conio.h>
#include<fstream>
#include<string>

using namespace std;

int main(){

	string fgicolor;
	string clcolor;
	cout << "--- PASTE THE CODE FROM trackList.txt IN YOUR tracklist.json FILE ---" << endl << endl;
	cout << "What is the new color of the fGI insert? (include # )" << endl;
	cin >> fgicolor;
	cout << "What is the new color of the CL? (include # )" << endl;
	cin >> clcolor;

	if (fgicolor.length() > 7 || clcolor.length() > 7){
		cout << "Invalid color entered." << endl;
		_getch();
		return 0;
	}
		
	ofstream file("consensusTrackList.txt");
	file << "{" << endl;
	file << " \"style\" : {  " << endl;
	file << " \t \"className\" : \"fGI_INS, CL\", " << endl;
	file << " \t \"color\" : \"function(feature, variableName, glyphObject, track) { if (feature.get('Type') == 'CL') { return '" << clcolor << "'; } else if (feature.get('Type') == 'fGI_INS') { return '" << fgicolor << "'; }  } \" , " << endl;
	file << " \t \"borderColor\" : \"function( feature, variableName, glyphObject, track ) { {return '#000000'; } } \" " << endl;
	file << " }, " << endl;
	file << " \"key\" : \"Consensus Track\", " << endl;
	file << " \"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", " << endl;
	file << " \"displayMode\" : \"collapsed\", " << endl;
	file << " \"trackType\" : \"null\", " << endl;
	file << " \"urlTemplate\" : \"tracks/Consensus Track/{refseq}/trackData.json\", " << endl;
	file << " \"compress\" : \"0\", " << endl;
	file << " \"label\" : \"Consensus Track\", " << endl;
	file << " \"type\" : \"CanvasFeatures\", " << endl;
	file << " \"category\" : \"Consensus\" " << endl;
	file << "}," << endl;

	cout << "Complete!  Check your consensusTrackList.txt file for the results. " << endl;
	_getch();
	return 0;


}