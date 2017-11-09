const int MAXV=300;

int cost[MAXV][MAXV], lvisit[MAXV], rvisit[MAXV];
int lvalue[MAXV], rvalue[MAXV], rmatch[MAXV], n;

inline bool findmatch(int s) {
	lvisit[s] = true;
	FOR(i,n) {
		if(rvisit[i]||lvalue[s]+rvalue[i]!=cost[s][i]) continue;
		rvisit[i] = true;
		if(rmatch[i]==-1 || findmatch(rmatch[i])) {
			rmatch[i] = s;
			return true;
		}
	}
	return false;
}

inline void normalize() {
	int down = INF;
	FOR(i,n) {
		if(lvisit[i]) {
			FOR(j,n) {
				if(!rvisit[j]) down = min(down,lvalue[i]+rvalue[j]-cost[i][j]);
			}
		}
	}
	FOR(i,n) {
		if(lvisit[i]) lvalue[i] -= down;
		if(rvisit[i]) rvalue[i] += down;
	}
}

inline int calc() {
	FOR(i,n) {
		lvalue[i] = *max_element(cost[i],cost[i]+n);
		rvalue[i] = 0;
	}
	FOR(i,n) rmatch[i] = -1;
	FOR(i,n) {
		memset(lvisit,0,sizeof lvisit);
		memset(rvisit,0,sizeof rvisit);
		if(findmatch(i)) continue;
		i--;
		normalize();
	}
	int ret = 0;
	FOR(i,n) ret += cost[rmatch[i]][i];
	return ret;
}
