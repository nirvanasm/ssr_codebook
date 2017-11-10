const int MAX_V=100;
struct Matching{
    int V,edge[MAX_V][MAX_V],match[MAX_V],dist[MAX_V],in_stack[MAX_V];
    vector<int> sta;
    Matching(int _V):V(_V){}
    void add_edge(int u,int v,int w){
        edge[u][v]=edge[v][u]=w;
    }
    bool SPFA(int u){
        if(in_stack[u])return 1;
        sta.push_back(u);
        in_stack[u]=1;
        for(int v=0;v<V;v++){
            if(u!=v&&match[u]!=v&&!in_stack[v]){
                int m=match[v];
                if(dist[m]>dist[u]-edge[v][m]+edge[u][v]){
                    dist[m]=dist[u]-edge[v][m]+edge[u][v];
                    sta.push_back(v);
                    in_stack[v]=1;
                    if(SPFA(m))return 1;
                    sta.pop_back();
                    in_stack[v]=0;
                }
            }
        }
        sta.pop_back();
        in_stack[u]=0;
        return 0;
    }
    int solve(){
        for(int i=0;i<V;i+=2){
            match[i]=i+1;
            match[i+1]=i;
        }
        while(1){
            bool found=0;
            for(int i=0;i<V;i++){
                dist[i]=in_stack[i]=0;
            }
            for(int i=0;i<V;i++){
                sta.clear();
                if(!in_stack[i]&&SPFA(i)){
                    found=1;
                    while(sta.size()>=2){
                        int u=sta.back();
                        sta.pop_back();
                        int v=sta.back();
                        sta.pop_back();
                        match[u]=v;
                        match[v]=u;
                    }
                }
            }
            if(!found)break;
        }
        int ret=0;
        for(int i=0;i<V;i++){
            ret+=edge[i][match[i]];
        }
        return ret/2;
    }
};
