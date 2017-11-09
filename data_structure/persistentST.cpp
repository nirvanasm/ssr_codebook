#include<bits/stdc++.h>
using namespace std;
struct nd{
	int val;
	nd *lf,*rg;
	nd(int x=0, nd* l=NULL,nd* r=NULL){val = x;lf = l;rg = r;}
};
nd* root[100005];
//BUILD ROOT[0]
nd* build(int l,int r){
	if(l==r)return new nd(0);
	int m = (l+r)/2;
	return new nd(0,build(l,m),build(m+1,r));
}
//ADD NUMBER X TO THE ARRAY
nd* upd(nd* bef, int l, int r, int x){
	if(l==r)return new nd(bef->val + 1);
	int m = (l+r)/2;
	if(x <= m)return new nd(bef->val + 1, upd(bef->lf,l,m,x), bef->rg);
	return new nd(bef->val + 1, bef->lf, upd(bef->rg,m+1,r,x));
}
//QUERY HOW MANY NUMBER LESS THAN MX AFTER N-TH UPDATE
int get(nd* cur, int l, int r, int mx){
	if(r<=mx)return cur->val;
	if(l>mx)return 0;	
	int m = (l+r)/2;
	return get(cur->lf, l, m, mx) + get(cur->rg, m+1, r, mx);
}
int main(){
	int n,m,ar[100005];//INPUT
	root[0]=build(0,n+5);
	for(int i = 0;i<n;i++)root[i+1] = upd(root[i],0,n+5,ar[i]);
}
