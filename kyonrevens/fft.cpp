typedef complex<double> Complex;
const double PI=acos(-1.0);
void fft(valarray<Complex> &v){
    int size=v.size();
    if(size<=1)return;
    valarray<Complex> even=v[slice(0,size/2,2)];
    valarray<Complex> odd=v[slice(1,size/2,2)];
    fft(even);
    fft(odd);
    for(int i=0;i<size/2;i++){
        Complex x=polar(1.0,-2*PI*i/size)*odd[i];
        v[i]=even[i]+x;
        v[i+size/2]=even[i]-x;
    }
}
void ifft(valarray<Complex> &v){
    v=v.apply(conj);
    fft(v);
    v=v.apply(conj);
    v/=v.size();
}
