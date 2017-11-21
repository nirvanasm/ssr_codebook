struct dt{
	int fi,se,idx;
	dt(int a=0,int b=0,int c=0){
		fi=a;se=b;idx=c;
	}
	bool operator != (dt x){
		return x.fi != this->fi || x.se != this->se;
	}
};

int ar[1005];
int sparse[LOG_MAX+5][1005];
string s;

void radix(vector<dt>&ls){
	vector<dt>tmp;tmp.resize(ls.size());
	int sz = ls.size();
	
	//SORT SECOND
	memset(ar,0,sizeof ar);
	for(int i = 0;i<sz;i++){
		ar[ls[i].se]++;
	}
	
	for(int i = 1;i<=sz;i++)ar[i]+=ar[i-1];
	for(int i = sz-1;i>=0;i--){
		tmp[--ar[ls[i].se]]=ls[i];
	}
	
	//SORT FIRST
	memset(ar,0,sizeof ar);
	for(int i = 0;i<sz;i++){
		ar[tmp[i].fi]++;
	}
	
	for(int i = 1;i<=sz;i++)ar[i]+=ar[i-1];
	for(int i = sz-1;i>=0;i--){
		ls[--ar[tmp[i].fi]]=tmp[i];
	}
}

int cnt[255];
void init(){
	memset(cnt,0,sizeof cnt);
	int mx = -1;
	for(int i=  0;i<s.length();i++){
		cnt[s[i]]++;
		mx = max(mx,(int)s[i]);
	}
	int cur = 1;
	for(int i = 0;i<=mx;i++){
		if(cnt[i]!=0)cnt[i]=cur++;
	}
	
	for(int i = 0;i<s.length();i++){
		sparse[0][i] = cnt[s[i]];
	}
}

int lcp(int id1,int id2){
	int ret = 0;
	for(int i = LOG_MAX, nx = (1<<10);i>=0;i--,nx/=2){
		if(id1 >= s.length() || id2 >= s.length())break;
		if(sparse[i][id1]==sparse[i][id2]){
			id1 += nx, id2 += nx;
			ret+=nx;
		}
	}
	return ret;
}

int main(){
  cin>>s;
  init();
  for(int i = 1,nx = 1;i<LOG_MAX;i++,nx*=2){
    vector<dt>ls;
    for(int j = 0;j<s.length();j++){
      if(j+nx<s.length())ls.pb(dt(sparse[i-1][j],sparse[i-1][j+nx],j));
      else ls.pb(dt(sparse[i-1][j],0,j));
    }
    radix(ls);
    int cur = 1;
    sparse[i][ls[0].idx]=cur;
    for(int j = 1;j<ls.size();j++){
      if(ls[j]!=ls[j-1])cur++;
      sparse[i][ls[j].idx]=cur;
    }
  }
}
