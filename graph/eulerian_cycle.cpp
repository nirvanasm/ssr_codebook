void fleurys(int idx){
  for(int i = 0;i<adj[idx].size();i++){
    int nex = adj[idx][i];
    if(mat[idx][nex]==0)continue;
    rmvEdge(idx,nex);
    dfs(nex);
  }
  ans.pb(idx);
}
//To find eulerian path simply make a dummy node for the starting and ending vertex
//i.e start DFS from start point or end point
