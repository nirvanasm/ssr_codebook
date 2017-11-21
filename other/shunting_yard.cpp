#include<bits/stdc++.h>
using namespace std;

const int BASE=10;

bool between(int x,int l,int r) {
	return (l<=x && x<=r);
}

stack<int> val;
stack<char> op;

inline void calculate() {
	int a=val.top(); val.pop();
	int b=val.top(); val.pop();
	char o=op.top(); op.pop();
	if (o=='+') val.push(a+b);
	if (o=='-') val.push(a-b);
	if (o=='*') val.push(a*b);
	if (o=='/') val.push(b/a);
}

bool precedence(char a,char b) {
	if (b=='('||b==')') return false;
	if ((a=='*'||a=='/')&&(b=='+'||b=='-')) return false;
	return true;
}

inline int evaluate(string &S) {
	for (int i=0;i<S.size();i++) {
		if (between(S[i],'0','9')) {
			int temp=0;
			while (i<S.size() && between(S[i],'0','9')) temp*=BASE, temp+=S[i++]-'0';
			i--; val.push(temp);
		}
		else if (S[i]=='(') op.push(S[i]);
		else if (S[i]==')') {
			while (op.top()!='(') calculate();
			op.pop();
		}
		else {
			while (!op.empty() && precedence(S[i],op.top())) calculate();
			op.push(S[i]);
		}
	}
	while (!op.empty()) calculate();
	int ret=val.top(); val.pop(); return ret;
}
