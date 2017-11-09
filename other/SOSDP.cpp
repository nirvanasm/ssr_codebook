FOR(mask,dua[n]){
	DP[mask]=A[0];
    for(int i=mask;i>0;i=(i-1)&mask) DP[mask]+=A[i];
}

FOR(i,dua[n]) DP[i]=A[i];
FOR(i,n) FOR(mask,dua[n]) if(mask&dua[i]) DP[mask]+=F[mask^dua[i]];
