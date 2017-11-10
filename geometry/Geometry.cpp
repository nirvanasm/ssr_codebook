/*General Double Operation*/

const double EPS=0.0;
const double PI=acos(-1.0);
const double INFD=1E9;

double between_d(double x,double l,double r) {
	return (min(l,r)<=x+EPS && x<=max(l,r)+EPS);
}

double same_d(double x,double y) {
	return between_d(x,y,y);
}

double dabs(double x) {
	if (x<EPS) return -x; return x;
}

int sign(double x) {
	return (0.0<x)-(x<0.0);
}

/*Point*/

struct point {
	double x,y;
	point() {
		x=y=0.0;
	}
	point(double _x,double _y) {
		x=_x; y=_y;
	}
	bool operator< (point other) const {
		if (y<other.y+EPS) return true;
		if (y+EPS>other.y) return false;
		return x<other.x+EPS;
	}
	bool operator== (point other) const {
		return same_d(x,other.x)&&same_d(y,other.y);
	}
};

double e_dist(point P1,point P2) {
	return hypot(P1.x-P2.x,P1.y-P2.y);
}

double m_dist(point P1,point P2) {
	return dabs(P1.x-P2.x)+dabs(P1.y-P2.y);
}

bool pointBetween(point P,point L,point R) {
	return (e_dist(L,P)+e_dist(P,R)==e_dist(L,R));
}

point mid(point P,point Q) {
	return point((P.x+Q.x)/2,(P.y+Q.y)/2);
}

/*Vector*/

struct vec {
	double x,y;
	vec() {
		x=y=0.0;
	}
	vec(double _x,double _y) {
		x=_x; y=_y;
	}
	vec(point A) {
		x=A.x; y=A.y;
	}
	vec(point A,point B) {
		x=B.x-A.x; y=B.y-A.y;
	}
};

vec scale(vec v,double s) {
	return vec(v.x*s,v.y*s);
}

vec flip(vec v) {
	return vec(-v.x,-v.y);
}

double dot(vec u,vec v) {
	return (u.x*v.x+u.y*v.y);
}

double cross(vec u,vec v) {
	return (u.x*v.y-u.y*v.x);
}

double norm_sq(vec v) {
	return (v.x*v.x+v.y*v.y);
}

point translate(point P,vec v) {
	return point(P.x+v.x,P.y+v.y);
}

point rotate(point P,point O,double angle) {
	vec v(O); P=translate(P,flip(v));
	return translate(point(P.x*cos(angle)-P.y*sin(angle),P.x*sin(angle)+P.y*cos(angle)),v);
}

double angle(point A,point O,point B) {
	vec OA(O,A), OB(O,B);
	return acos(dot(OA,OB)/sqrt(norm_sq(OA)*norm_sq(OB)));
}

double orientation(point O,point P,point Q) {
	vec OP(O,P), OQ(O,Q);
	double c=cross(OP,OQ);
	return c;
}

/*Line*/

struct line {
	double a,b,c;
	line() {
		a=b=c=0.0;
	}
	line(double _a,double _b,double _c) {
		a=_a; b=_b; c=_c;
	}
	line(point P1,point P2) {
		if (P2<P1) {
			point T; T=P1; P1=P2; P2=T;
		}
		if (same_d(P1.x,P2.x)) a=1.0, b=0.0, c=-P1.x;
		else a=-(P1.y-P2.y)/(P1.x-P2.x), b=1.0, c=-(a*P1.x)-P1.y;
	}
	line (point P,double slope) {
		if (same_d(slope,INFD)) a=1.0, b=0.0, c=-P.x;
		else a=-slope, b=1.0, c=-(a*P.x)-P.y;
	}
	bool operator== (line other) const {
		return same_d(a,other.a)&&same_d(b,other.b)&&same_d(c,other.c);
	}
	double slope() {
		if (same_d(b,0.0)) return INFD;
		return -(a/b);
	}
};

bool paralel(line L1,line L2) {
	return same_d(L1.a,L2.a)&&same_d(L1.b,L2.b);
}

bool intersection(line L1,line L2,point &P) {
	if (paralel(L1,L2)) return false;
	P.x=(L2.b*L1.c-L1.b*L2.c)/(L2.a*L1.b-L1.a*L2.b);
	if (same_d(L1.b,0.0)) P.y=-(L2.a*P.x+L2.c);
	else P.y=-(L1.a*P.x+L1.c);
	return true;
}

double pointToLine(point P,point A,point B,point &C) {
	vec AP(A,P), AB(A,B);
	double u=dot(AP,AB)/norm_sq(AB);
	C=translate(A,scale(AB,u));
	return e_dist(P,C);
}

double lineToLine(line L1,line L2) {
	if (!paralel(L1,L2)) return 0.0;
	return dabs(L2.c-L1.c)/sqrt(L1.a*L1.a+L1.b*L1.b);
}

/*Line Segment*/

struct segment {
	point P,Q;
	line L;
	segment() {
		point T1; P=Q=T1;
		line T2; L=T2;
	}
	segment(point _P,point _Q) {
		if (_Q<_P) {
			point T1=_P; _P=_Q; _Q=T1;
		}
		P=_P; Q=_Q;
		line T2(_P,_Q); L=T2;
	}
	bool operator== (segment other) const {
		return P==other.P&&Q==other.Q;
	}
};

bool onSegment(point P,segment S) {
	if (orientation(S.P,S.Q,P)!=0.0) return false;
	return between_d(P.x,S.P.x,S.Q.x) && between_d(P.y,S.P.y,S.Q.y);
}

bool s_intersection(segment S1,segment S2) {
	double o1=orientation(S1.P,S1.Q,S2.P);
	double o2=orientation(S1.P,S1.Q,S2.Q);
	double o3=orientation(S2.P,S2.Q,S1.P);
	double o4=orientation(S2.P,S2.Q,S1.Q);
	if (o1!=o2 && o3!=o4) return true;
	if (o1==0.0 && onSegment(S2.P,S1)) return true;
	if (o2==0.0 && onSegment(S2.Q,S1)) return true;
	if (o3==0.0 && onSegment(S1.P,S2)) return true;
	if (o4==0.0 && onSegment(S2.Q,S2)) return true;
	return false;
}

int ss_intersection(segment S1,segment S2,point &P1,point &P2) {
	if (S2.P<S1.P) return ss_intersection(S2,S1,P1,P2);
	if (intersection(S1.L,S2.L,P1)) {
		if (onSegment(P1,S1) && onSegment(P1,S2)) return 1;
		return 0;
	}
	if (!(S1.L==S2.L)) return 0;
	if (S1.Q==S2.P) {
		P1=S1.Q; return 1;
	}
	P1=S2.P; P2=S1.Q; return 2;
}

double pointToSegment(point P,point A,point B,point &C) {
	vec AP(A,P), AB(A,B);
	double u=dot(AP,AB)/norm_sq(AB);
	if (u<EPS) {
		C=A; return e_dist(P,A);
	}
	if (u+EPS>1.0) {
		C=B; return e_dist(P,B);
	}
	return pointToLine(P,A,B,C);
}

double segmentToSegment(segment S1,segment S2) {
	if (s_intersection(S1,S2)) return 0.0;
	double ret=INFD; point dummy;
	ret=min(ret,pointToSegment(S1.P,S2.P,S2.Q,dummy));
	ret=min(ret,pointToSegment(S1.Q,S2.P,S2.Q,dummy));
	ret=min(ret,pointToSegment(S2.P,S1.P,S2.Q,dummy));
	ret=min(ret,pointToSegment(S2.Q,S1.P,S2.Q,dummy));
	return ret;
}

/*Circle*/

struct circle {
	point P;
	double r;
	circle() {
		point P1; P=P1;
		r=0.0;
	}
	circle(point _P,double _r) {
		P=_P; r=_r;
	}
	circle(point P1,point P2) {
		P=mid(P1,P2); r=e_dist(P,P1);
	}
	circle(point P1,point P2,point P3) {
		vector<point> T; T.clear(); T.pb(P1); T.pb(P2); T.pb(P3); sort(T.begin(),T.end());
		P1=T[0]; P2=T[1]; P3=T[2];
		point M1,M2; M1=mid(P1,P2); M2=mid(P2,P3);
		point Q2,Q3; Q2=rotate(P2,P1,PI/2); Q3=rotate(P3,P2,PI/2);
		vec P1Q2(P1,Q2), P2Q3(P2,Q3);
		point M3,M4; M3=translate(M1,P1Q2); M4=translate(M2,P2Q3);
		line L1(M1,M3), L2(M2,M4);
		intersection(L1,L2,P); r=e_dist(P,P1);
	}
	bool operator==(circle other) const {
		return (P==other.P && same_d(r,other.r));
	}
};

bool insideCircle(point P,circle C) {
	return e_dist(P,C.P)<=C.r+EPS;
}

bool c_intersection(circle C1,circle C2,point &P1,point &P2) {
	double d=e_dist(C1.P,C2.P);
	if (d+EPS>C1.r+C2.r) return false;
	if (d<dabs(C1.r-C2.r)+EPS) return false;
	double x1=C1.P.x, y1=C1.P.y, r1=C1.r, x2=C2.P.x, y2=C2.P.y, r2=C2.r;
	double a=(r1*r1-r2*r2+d*d)/(2*d), h=sqrt(r1*r1-a*a);
	point T(x1+a*(x2-x1)/d,y1+a*(y2-y1)/d);
	P1=point(T.x-h*(y2-y1)/d,T.y+h*(x2-x1)/d); P2=point(T.x+h*(y2-y1)/d,T.y-h*(x2-x1)/d);
	return true;
}

bool lc_intersection(line L,circle O,point &P1,point &P2) {
	double a=L.a, b=L.b, c=L.c, x=O.P.x, y=O.P.y, r=O.r;
	double A=a*a+b*b, B=2*a*b*y-2*a*c-2*b*b*x, C=b*b*x*x+b*b*y*y-2*b*c*y+c*c-b*b*r*r;
	double D=B*B-4*A*C; point T1,T2;
	if (same_d(b,0.0)) {
		T1.x=c/a;
		if (dabs(x-T1.x)+EPS>r) return false;
		if (same_d(T1.x-r-x,0.0) || same_d(T1.x+r-x,0.0)) {
			P1=P2=point(T1.x,y); return true;
		}
		double dx=dabs(T1.x-x), dy=sqrt(r*r-dx*dx);
		P1=point(T1.x,y-dy); P2=point(T1.x,y+dy); return true;
	}
	if (same_d(D,0.0)) {
		T1.x=-B/(2*A); T1.y=(c-a*T1.x)/b; P1=P2=T1; return true;
	}
	if (D<EPS) return false;
	D=sqrt(D);
	T1.x=(-B-D)/(2*A);T1.y=(c-a*T1.x)/b; P1=T1;
	T2.x=(-B+D)/(2*A);T2.y=(c-a*T2.x)/b; P2=T2; return true;
}

bool sc_intersection(segment S,circle C,point &P1,point &P2) {
	bool cek=lc_intersection(S.L,C,P1,P2);
	if (!cek) return false;
	double x1=S.P.x, y1=S.P.y, x2=S.Q.x, y2=S.Q.y;
	bool b1=between_d(P1.x,x1,x2)&&between_d(P1.y,y1,y2);
	bool b2=between_d(P2.x,x1,x2)&&between_d(P2.y,y1,y2);
	if (P1==P2) return b1;
	if (b1||b2) {
		if (!b1) P1=P2; if (!b2) P2=P1; return true;
	}
	return false;
}

/*Triangle*/

double t_perimeter(point A,point B,point C) {
	return e_dist(A,B)+e_dist(B,C)+e_dist(C,A);
}

double t_area(point A,point B,point C) {
	double s=t_perimeter(A,B,C)/2;
	double ab=e_dist(A,B), bc=e_dist(B,C), ac=e_dist(C,A);
	return sqrt(s*(s-ab)*(s-bc)*(s-ac));
}

circle t_inCircle(point A,point B,point C) {
	vector<point> T; T.clear(); T.pb(A); T.pb(B); T.pb(C); sort(T.begin(),T.end());
	A=T[0]; B=T[1]; C=T[2];
	double r=t_area(A,B,C)/(t_perimeter(A,B,C)/2);
	double ratio=e_dist(A,B)/e_dist(A,C);
	vec BC(B,C); BC=scale(BC,ratio/(1+ratio));
	point P; P=translate(B,BC); line AP1(A,P);
	ratio=e_dist(B,A)/e_dist(B,C);
	vec AC(A,C); AC=scale(AC,ratio/(1+ratio));
	P=translate(A,AC); line BP2(B,P);
	intersection(AP1,BP2,P); return circle(P,r);
}

circle t_outCircle(point A,point B,point C) {
	return circle(A,B,C);
}

/*Polygon CW*/

struct polygon {
	vector<point> P;
	polygon() {
		P.clear();
	}
	polygon(vector<point> &_P) {
		P=_P;
	}
	int prev(int idx) {
		return (idx==0?P.size()-1:idx-1);
	}
	int next(int idx) {
		return (idx==P.size()-1?0:idx+1);
	}
	double perimeter() {
		double ret=0;
		FOR(i,P.size()) {
			ret+=e_dist(P[i],P[next(i)]);
		}
		return ret;
	}
	double area() {
		double ret=0;
		FOR(i,P.size()) {
			ret+=P[i].x*(P[prev(i)].y-P[next(i)].y);
		}
		return ret/2;
	}
};

polygon convexHull(vector<point> &pts){
    sort(pts.begin(),pts.end());
    vector<point> hull;
    for(int i=0;i<2;i++){
        int start=(int)hull.size();
        for(auto pt:pts){
            while((int)hull.size()>=start+2&&orientation(hull[(int)hull.size()-1],hull[(int)hull.size()-2],pt)<=0.0)hull.pob();
            hull.pb(pt);
        }
        hull.pop_back();
        reverse(pts.begin(),pts.end());
    }
    if((int)hull.size()==2&&hull[0]==hull[1])hull.pob();
    return polygon(hull);
}

int findTop(polygon &A) {
	point P=A.P[0];
	int ret=0;
	FOR(i,A.P.size()) {
		if (P<A.P[i]) P=A.P[i], ret=i;
	}
	return ret;
}

int inConvexPolygon(point &P,polygon &A,int idxTop) { //-1 inside,0 onsegment, 1 outside
	if (P<A.P[0]||A.P[idxTop]<P) return 1;
	double o=orientation(P,A.P[0],A.P[idxTop]);
	if (o==0.0) {
		if (P==A.P[0]||P==A.P[idxTop]) return 0;
		return ((idxTop==1||idxTop+1==(int)A.P.size())?0:-1);
	}
	else if (o<0.0) {
		vector<point>::reverse_iterator itLeft=upper_bound(A.P.rbegin(),A.P.rend()-idxTop-1,P);
		return sign(orientation((itLeft==A.P.rbegin())?A.P[0]:itLeft[-1],P,itLeft[0]));
	}
	else {
		vector<point>::iterator itRight=lower_bound(A.P.begin()+1,A.P.begin()+idxTop,P);
		return sign(orientation(itRight[0],P,itRight[-1]));
	}
}

int inSimplePolygon(point &P,polygon &A) { //-1 inside,0 onsegment, 1 outside
	int ret=0;
	FOR(i,A.P.size()) {
		if (P==A.P[i]) return 0;
		int j=A.next(i);
		if (onSegment(P,segment(A.P[i],A.P[j]))) return 0;
		bool below=(A.P[i].y<P.y);
		if (below!=(A.P[j].y<P.y)) {
			double o=orientation(P,A.P[i],A.P[j]);
			if (o==0.0) return 0;
			if (below==(o>0.0)) ret+=below?1:-1;
		}
	}
	return ret==0?1:-1;
}

circle minCoverCircle(polygon &A) {
	vector<point> p=A.P;
    point c; circle ret;
    double cr = 0.0;
    int i, j, k;
    c = p[0];
    for(i = 1; i < p.size(); i++) {
        if(e_dist(p[i], c) >= cr+EPS) {
            c = p[i], cr = 0;
            for(j = 0; j < i; j++) {
                if(e_dist(p[j], c) >= cr+EPS) {
                	c=mid(p[i],p[j]);
                    cr = e_dist(p[i], c);
                    for(k = 0; k < j; k++) {
                        if(e_dist(p[k], c) >= cr+EPS) {
                        	ret=circle(p[i],p[j],p[k]);
                            c=ret.P; cr=ret.r;
                        }
                    }
                }
            }
        }
    }
    return ret;
}

/*Geometry Algorithm*/

double DP[110][110];
double minCostPolygonTriangulation(polygon &A) {
	if (A.P.size()<3) return 0;
	REP(i,0,A.P.size()) {
   		for (int j=0,k=i;k<A.P.size();j++,k++) {
   			if (k<j+2) DP[j][k]=0.0;
          	else {
            	DP[j][k]=INFD;
            	REP(l,j+1,k) {
            		double cost=e_dist(A.P[j],A.P[k])+e_dist(A.P[k],A.P[l])+e_dist(A.P[l],A.P[j]);
            		DP[j][k]=min(DP[j][k],DP[j][l]+DP[l][k]+cost);
				}
			}
		}
    }
	return DP[0][A.P.size()-1];
}
