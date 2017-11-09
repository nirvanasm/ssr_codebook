int f[1000010];

inline void precomp(string S) {
	int len=0,i=1; f[0]=0;
	while (i<S.size()) {
		if (S[i]==S[len]) f[i++]=++len;
		else if (len!=0) len=f[len-1];
		else f[i++]=0;
	}
}

vector<int> search(string txt,string pat) {
	vector<int> ret; ret.clear();
	precomp(pat);
	int idxT,idxP; idxT=idxP=0;
	while (idxT<txt.size()) {
		if (txt[idxT]==pat[idxP]) idxT++, idxP++;
		if (idxP==pat.size()) ret.pb(idxT-idxP), idxP=f[idxP-1];
		else if (idxT<txt.size() && txt[idxT]!=pat[idxP]) {
			if (idxP!=0) idxP=f[idxP-1];
			else idxT++;
		}
	}
	return ret;
}
