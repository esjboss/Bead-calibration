%function bhmie calculates amplitude scatteing matrix elements and efficiencies for extinction, 
%total scattering and backscattering for a given size parameter and relative refractive index
%this code is a translation of a Fortran code of Bohren & Huffman
%Emmanuel Boss Sep 10 1998

function[S1,S2,Qscat,Qext,Qback]=bhmie(x,refrel,nang);

dx=x;
y=x*refrel;

%series terminated after nstop terms

nstop=ceil(x+4*x^.3333+2);
ymod=abs(y);
nmx=max(nstop,ceil(ymod))+15;
dang=pi/2/(nang-1);
for j=1:nang
   theta(j)=(j-1)*dang;
   AMU(j)=cos(theta(j));
end
%Logarithmic derivative D(j) calculated by downward recurence 
   
D(nmx)=0+0i;
nn=nmx-1;
for n=1:nn
   RN=nmx-n+1;
   D(nmx-n)=(RN/y)-(1./(D(nmx-n+1)+RN/y));
end
PI0=zeros(nang,1);
PI1=ones(nang,1);
nn=2*nang-1;
S1=zeros(nn,1);
S2=zeros(nn,1);
%Riccati-Bessel functins with real argument x calculated by upward recurence
PSI0=cos(dx);
PSI1=sin(dx);
CHI0=-sin(x);
CHI1=cos(x);
APSI0=PSI0;
APSI1=PSI1;
XI0=APSI0-CHI0*sqrt(-1);
XI1=APSI1-CHI1*sqrt(-1);
Qscat=0.0;
n=1;
while(n-1-nstop<0)
	DN=n;
   RN=n;
   FN=(2*RN+1)/(RN*(RN+1));
   PSI=(2*DN-1)*PSI1/dx-PSI0;
   APSI=PSI;
   CHI=(2*RN-1)*CHI1/x-CHI0;
   XI=APSI-CHI*sqrt(-1);
   AN=((D(n)/refrel+RN/x)*APSI-APSI1)/((D(n)/refrel+RN/x)*XI-XI1);
   BN=((refrel*D(n)+RN/x)*APSI-APSI1)/((D(n)*refrel+RN/x)*XI-XI1);
   Qscat=Qscat+(2.*RN+1)*(abs(AN)^2+abs(BN)^2);
   for j=1:nang
      jj=2*nang-j;
      PI(j)=PI1(j);
      TAU(j)=RN*AMU(j)*PI(j)-(RN+1)*PI0(j);
      P=(-1)^(n-1);
      S1(j)=S1(j)+FN*(AN*PI(j)+BN*TAU(j));
      T=(-1)^n;
      S2(j)=S2(j)+FN*(AN*TAU(j)+BN*PI(j));
      if(j~=jj)
     	S1(jj)=S1(jj)+FN*(AN*PI(j)*P+BN*TAU(j)*T);
      	S2(jj)=S2(jj)+FN*(BN*PI(j)*P+AN*TAU(j)*T);
	  end
   end
   PSI0=PSI1;
   PSI1=PSI;
   APSI1=PSI1;
   CHI0=CHI1;
   CHI1=CHI;
   XI1=APSI1-CHI1*sqrt(-1);
   n=n+1;
   RN=n;
   for j=1:nang
      PI1(j)=((2*RN-1)/(RN-1))*AMU(j)*PI(j)-RN*PI0(j)/(RN-1);
      PI0(j)=PI(j);
   end
end
Qscat=(2./x^2)*Qscat;
Qext=(4./x^2)*real(S1(1));
%Note: Qback is not Qbb, but the radar back scattering.
Qback=(4./x^2)*abs(S1(2*nang-1))^2;
