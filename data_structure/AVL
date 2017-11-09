struct node {
	int key;
	int height;
	node* left;
	node* right;
};

int getHeight(node* x) {
	if (x==NULL) return 0;
	return x->height;
}

int getBalance(node* x) {
	if (x==NULL) return 0;
	return getHeight(x->left)-getHeight(x->right);
}

void updateHeight(node* x) {
	x->height=max(getHeight(x->left),getHeight(x->right))+1;
}

node* newNode(int key) {
	node* ret=new node;
	ret->key=key; ret->height=1;
	ret->left=ret->right=NULL;
	return ret;
}

node* rightRotate(node* y) {
	node* x=y->left;
	node* B=x->right;
	x->right=y; y->left=B;
	updateHeight(y); updateHeight(x); return x;
}

node* leftRotate(node* x) {
	node* y=x->right;
	node* B=y->left;
	y->left=x; x->right=B;
	updateHeight(x); updateHeight(y); return y;
}

node* insert(node* cur,int key) {
	if (cur==NULL) return newNode(key);
	if (key==cur->key) return cur;
	if (key<cur->key) cur->left=insert(cur->left,key);
	else cur->right=insert(cur->right,key);
	updateHeight(cur);
	int balance=getBalance(cur);
	if (balance>1) {
		if (key>cur->left->key) cur->left=leftRotate(cur->left);
		return rightRotate(cur);
	}
	if (balance<-1) {
		if (key<cur->right->key) cur->right=rightRotate(cur->right);
		return leftRotate(cur);
	}
	return cur;
}

node* findData(node* cur,int key) {
	if (cur==NULL) return NULL;
	if (key==cur->key) return cur;
	if (key<cur->key) return findData(cur->left,key);
	return findData(cur->right,key);
}

node* successor(node* x) {
	node* cur=x->right;
	while (cur->left!=NULL) cur=cur->left;
	return cur;
}

node* clear(node* cur) {
	if (cur!=NULL) {
		clear(cur->left); clear(cur->right); delete cur;
	}
	return NULL;
}

node* deleteValue(node* cur,int key) {
	if (cur==NULL) return NULL;
	if (key<cur->key) cur->left=deleteValue(cur->left,key);
	else if (key>cur->key) cur->right=deleteValue(cur->right,key);
	else {
		if (cur->left==NULL || cur->right==NULL) {
			node* temp=NULL;
			if (cur->left!=NULL) temp=cur->left;
			else if (cur->right!=NULL) temp=cur->right;
			if (temp==NULL) temp=cur, cur=NULL;
			else *cur=*temp;
			delete temp;
		}
		else {
			node* temp=successor(cur);
			cur->key=temp->key;
			cur->right=deleteValue(cur->right,temp->key);
		}
	}
	if (cur==NULL) return NULL;
	updateHeight(cur);
	int balance=getBalance(cur);
	if (balance>1 && getBalance(cur->left)>=0) return rightRotate(cur);
	if (balance>1 && getBalance(cur->left)<0) {
		cur->left=leftRotate(cur->left); return rightRotate(cur);
	}
	if (balance<-1 && getBalance(cur->right)<=0) return leftRotate(cur);
	if (balance<-1 && getBalance(cur->right)>0) {
		cur->right=rightRotate(cur->right); return leftRotate(cur);
	}
	return cur;
}

void print(node* cur) {
	if (cur==NULL) return;
	print(cur->left);
	cout<<"> "<<cur->key<<endl;
	print(cur->right);
}

node* AVL=NULL;
