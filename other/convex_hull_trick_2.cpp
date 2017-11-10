struct Line{
    LL M,C;
    Line(){}
    Line(LL _M,LL _C):M(_M),C(_C){}
};
struct ConvexHull{
    int pointer;
    vector<Line> L;
    ConvexHull():pointer(0){}
    bool check(Line a,Line b,Line c){
        return (a.C-b.C)*(c.M-a.M)>(a.C-c.C)*(b.M-a.M);
    }
    void add(Line l){
        while(L.size()>=2&&!check(L[L.size()-2],L[L.size()-1],l)){
            if(pointer==L.size()-1)pointer--;
            L.pop_back();
        }
        L.push_back(l);
    }
    LL query(LL x){
        while(pointer+1<L.size()&&L[pointer].M*x+L[pointer].C>L[pointer+1].M*x+L[pointer+1].C)pointer++;
        return L[pointer].M*x+L[pointer].C;
    }
};
