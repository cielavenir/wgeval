/*
wgeval: yet another evaluator of wgsim/DWGSIM
(C) T. Yamada under 2-clause BSDL.

usage:
wgeval < mapping.sam > mapping.wgeval
tail -n1 mapping.wgeval # shows the number of correctly mapped reads
wc -l mapping.wgeval # shows the number of wrongly mapped reads (subtract the trailing last line)
*/

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
using namespace std;

typedef map<string,bool> MAP;

vector<string> split(string &str, const char *delim){
	vector<string> result;
	int cutAt;
	while( (cutAt = str.find_first_of(delim)) != str.npos ){
		if(cutAt > 0){
			result.push_back(str.substr(0, cutAt));
		}
		str = str.substr(cutAt + 1);
	}
	if(str.length() > 0){
		result.push_back(str);
	}
	return result;
}

int main(int argc, char **argv){
	cin.tie(0);
	ios::sync_with_stdio(false);
	int gap=0;
	MAP m;
	string s;
	for(;getline(cin,s);){
		if(s[0]=='@')continue;
		int i;
		vector<string> t=split(s,"\t");
		if(t.size()<11)continue; //something is wrong.
		bool correct=true;

		//parameters
		int flag=strtol(t[1].c_str(),NULL,10);
		string chromo=t[2];
		long long left=strtoll(t[3].c_str(),NULL,10);
		long long right=strtoll(t[3].c_str(),NULL,10);
		int quality=strtol(t[4].c_str(),NULL,10)/10;
		if((flag&0x4)||chromo=="*")continue; //unmapped

		//cigar
		string cigar=t[5];
		char *rem=NULL;
		int length=strtol(cigar.c_str(),&rem,10);
		if(rem){
			for(i=0;i<strlen(rem);i++)if(rem[i]<'A'||'Z'<rem[i])break; //fixme
			if(i<strlen(rem)){cerr<<"Cannot parse cigar: "<<cigar<<endl;}//return 1;}

			for(i=0;i<strlen(rem);i++){
				if(rem[i]=='M'||rem[i]=='D'||rem[i]=='N')right+=length;
			}
		}
		right--;
		long long left0=left,right0=right;
		if(rem){
			for(i=0;i<strlen(rem);i++){
				if(rem[i]=='S'||rem[i]=='H'){
					left-=length;
					right+=length;
					left0-=length;
					right0+=length;
				}
			}
		}

		//read name
		string read_name=t[0];
		if(t[0][t[0].size()-2]=='/'){ //wgsim
			vector<string> read_args=split(t[0],"_");
			if(read_args[0]=="rand")continue;
			if(read_args[0]!=chromo)correct=false;
			else if(flag&0x10){ //reverse
				long long coor=strtoll(read_args[2].c_str(),NULL,10);
				if(abs(coor - right) > gap && abs(coor - right0) > gap)correct=false;
			}else{
				long long coor=strtoll(read_args[1].c_str(),NULL,10);
				if(abs(coor - left) > gap && abs(coor - left0) > gap)correct=false;
			}
		}else{ //dwgsim
			vector<string> read_args=split(t[0],"_");
			long long coor=strtoll(read_args[(flag&0x80)?2:1].c_str(),NULL,10);
			if(read_args[0]=="rand")continue;
			if(read_args[0]!=chromo)correct=false;
			if(abs(coor - left) > gap && abs(coor - left0) > gap)correct=false;
		}
		m[read_name]|=correct;
	}
	MAP::iterator it=m.begin();
	int n=0;
	for(;it!=m.end();it++){
		if(!it->second)cout<<it->first<<endl;
		else n++;
	}
	cout<<n<<endl;
}