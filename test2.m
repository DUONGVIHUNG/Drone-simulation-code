%Physical parameter
syms M R m l %kitai 
syms mi a b  %prop
syms x y z gamma beta alpha th1 th2 th3 th4
Ixx=@(M,R,m,l)(0.4*M*R^2+2*m*l^2);
Izz=@(M,R,m,l)(0.5*M*R^2+4*m*l^2);
Ix=Ixx(M,R,m,l);
Iy=Ix;
Iz=Izz(M,R,m,l);
Io=[Ix 0  0;...
    0  Iy 0;...
    0  0 Iz;];
Mo=[(M+m)*eye(3) zeros(3);...
    zeros(3)             Io];

Ixxi=@(m,b)(1/12*m*b^2);
Iyyi=@(m,a)(1/12*m*a^2);

Ix1=Ixxi(mi,b);
Iy1=Iyyi(mi,a);
Iz1=simplify(Ix1+Iy1);
I1=[Ix1 0  0;...
    0  Iy1 0;...
    0   0 Iz1;];
M1=[mi*eye(3)  zeros(3);...
    zeros(3)     I1];
M2=M1;
M3=M1;
M4=M1;
%% Exponential parameter
% imagine link
qx=[0;0;0];
qy=qx;
qz=qx;
qg=qx;%gamma
qb=qx;%beta
qa=qx;%alpha
% props as link
q1=[l*sqrt(2)/2;l*sqrt(2)/2;0];
q2=[l*sqrt(2)/2;-l*sqrt(2)/2;0];
q3=[-l*sqrt(2)/2;-l*sqrt(2)/2;0];
q4=[-l*sqrt(2)/2;l*sqrt(2)/2;0];
% translation joint parameter
vx=[1;0;0];
vy=[0;1;0];
vz=[0;0;1];
% rotaion joint parameter
wg=[1;0;0];
wb=[0;1;0];
wa=[0;0;1];
w1=wa;
w2=wa;
w3=wa;
w4=wa;

%% Twist
% for translation
twt=@(v)...
    ([v;zeros(3,1)]);
twx=twt(vx);
twy=twt(vy);
twz=twt(vz);
% for rotation
twr=@(w,q)...
    ([-cross(w,q);w]);
twg=twr(wg,qg);
twb=twr(wb,qb);
twa=twr(wa,qa);
tw1=twr(w1,q1);
tw2=twr(w2,q2);
tw3=twr(w3,q3);
tw4=twr(w4,q4);
%% calc for exponential of all twist
% for translation
et=@(v,th)...
    ([eye(3)   v*th;...
    zeros(1,3)  1;]);
ex=et(vx,x);
ey=et(vy,y);
ez=et(vz,z);
% for rotation
wh=@(w)... % transform vector to skew matrix
    ([0  -w(3,1)   w(2,1);...
      w(3,1) 0     w(1,1);...
     -w(2,1) -w(1,1)  0]);
whg=wh(wg);
whb=wh(wb);
wha=wh(wa);
wh1=wh(w1);
wh2=wh(w2);
wh3=wh(w3);
wh4=wh(w4);
ew=@(wh,th)...% calc for exponential of omega
    (eye(3)+sin(th)*wh+wh*wh*(1-cos(th)));
ewg=ew(whg,gamma);
ewb=ew(whb,beta);
ewa=ew(wha,alpha);
ew1=ew(wh1,th1);
ew2=ew(wh2,th2);
ew3=ew(wh3,th3);
ew4=ew(wh4,th4);
er=@(ew,w,q)... @calc for exponential of twist(e^guzai)
    ([eye(3)      (eye(3)-ew)*cross(w,-cross(w,q));...
     zeros(1,3)                1                  ]);
eg=er(ewg,wg,qg);
eb=er(ewb,wb,qb);
ea=er(ewa,wa,qa);
e1=er(ew1,w1,q1);
e2=er(ew2,w2,q2);
e3=er(ew3,w3,q3);
e4=er(ew4,w4,q4);
%% Solve for M,C,N
% calculate for link configuration when all theta=0
gsl=@(R,p)...
    ([R             p;...
   zeros(1,3)       1]);
% kitai config
gslo=gsl(eye(3),qx);
% props config
gsl1=gsl(eye(3),q1);
gsl2=gsl(eye(3),q2);
gsl3=gsl(eye(3),q3);
gsl4=gsl(eye(3),q4);
% production of exponential 
h=@(x,y,z,g,b,a,t,gs)...
    (x*y*z*g*b*a*t*gs);
% production of exponential for body(frame) link
% for frame
hox=h(ex,ey,ez,eg,eb,ea,1,gslo);
hoy=h(1,ey,ez,eg,eb,ea,1,gslo);
hoz=h(1,1,ez,eg,eb,ea,1,gslo);
hog=h(1,1,1,eg,eb,ea,1,gslo);
hob=h(1,1,1,1,eb,ea,1,gslo); 
hoa=h(1,1,1,1,1,ea,1,gslo);
% for prop 
%calc for production of exponential
poex=ex*ey*ez*eg*eb*ea;
poey=ey*ez*eg*eb*ea;
poez=ez*eg*eb*ea;
poeg=eg*eb*ea;
poeb=eb*ea;
poea=ea;
% calc for prop config
hi=@(poe,e,gs)(poe*e*gs);
  %for prop 1
  h1x=hi(poex,e1,gsl1);
  h1y=hi(poey,e1,gsl1);
  h1z=hi(poez,e1,gsl1);
  h1g=hi(poeg,e1,gsl1);
  h1b=hi(poeb,e1,gsl1);
  h1a=hi(poea,e1,gsl1);
  h1=hi(1,e1,gsl1);
  %for prop 2
  h2x=hi(poex,e2,gsl2);
  h2y=hi(poey,e2,gsl2);
  h2z=hi(poez,e2,gsl2);
  h2g=hi(poeg,e2,gsl2);
  h2b=hi(poeb,e2,gsl2);
  h2a=hi(poea,e2,gsl2);
  h2=hi(1,e2,gsl2);
  %for prop 3
  h3x=hi(poex,e3,gsl3);
  h3y=hi(poey,e3,gsl3);
  h3z=hi(poez,e3,gsl3);
  h3g=hi(poeg,e3,gsl3);
  h3b=hi(poeb,e3,gsl3);
  h3a=hi(poea,e3,gsl3);
  h3=hi(1,e3,gsl3);
  %for prop 4
  h4x=hi(poex,e4,gsl4);
  h4y=hi(poey,e4,gsl4);
  h4z=hi(poez,e4,gsl4);
  h4g=hi(poeg,e4,gsl4);
  h4b=hi(poeb,e4,gsl4);
  h4a=hi(poea,e4,gsl4);
  h4=hi(1,e4,gsl4);

% extract R from h
R=eye(3);
% extract p from h
p=@(h)...
    (h(1:3,4));
   %for frame
   pox=p(hox);
   poy=p(hoy);
   poz=p(hoz);
   pog=p(hog);
   pob=p(hob);
   poa=p(hoa);
   %for prop 1
   p1x=p(h1x);
   p1y=p(h1y);
   p1z=p(h1z);
   p1g=p(h1g);
   p1b=p(h1b);
   p1a=p(h1a);
   p1=p(h1);
   %for prop 2
   p2x=p(h2x);
   p2y=p(h2y);
   p2z=p(h2z);
   p2g=p(h2g);
   p2b=p(h2b);
   p2a=p(h2a);
   p2=p(h2);
   %for prop 3
   p3x=p(h3x);
   p3y=p(h3y);
   p3z=p(h3z);
   p3g=p(h3g);
   p3b=p(h3b);
   p3a=p(h3a);
   p3=p(h3);
   %for prop 4
   p4x=p(h4x);
   p4y=p(h4y);
   p4z=p(h4z);
   p4g=p(h4g);
   p4b=p(h4b);
   p4a=p(h4a);
   p4=p(h4);
% calc for inverse adjoint matrix of h
%function
Adi=@(R,p)...
    ([R.'   -R.'*wh(p);...
   zeros(3)        R.']);
%for frame
 Adix=Adi(R,pox);
 Adiy=Adi(R,poy);
 Adiz=Adi(R,poz);
 Adig=Adi(R,pog);
 Adib=Adi(R,pob);
 Adia=Adi(R,poa);
%for prop 1
Adi1x=Adi(R,p1x);
Adi1y=Adi(R,p1y);
Adi1z=Adi(R,p1z);
Adi1g=Adi(R,p1g);
Adi1b=Adi(R,p1b);
Adi1a=Adi(R,p1a);
Adi1=Adi(R,p1);
%for prop 2
Adi2x=Adi(R,p2x);
Adi2y=Adi(R,p2y);
Adi2z=Adi(R,p2z);
Adi2g=Adi(R,p2g);
Adi2b=Adi(R,p2b);
Adi2a=Adi(R,p2a);
Adi2=Adi(R,p2);
%for prop 3
Adi3x=Adi(R,p3x);
Adi3y=Adi(R,p3y);
Adi3z=Adi(R,p3z);
Adi3g=Adi(R,p3g);
Adi3b=Adi(R,p3b);
Adi3a=Adi(R,p3a);
Adi3=Adi(R,p3);
% for prop 4
Adi4x=Adi(R,p4x);
Adi4y=Adi(R,p4y);
Adi4z=Adi(R,p4z);
Adi4g=Adi(R,p4g);
Adi4b=Adi(R,p4b);
Adi4a=Adi(R,p4a);
Adi4=Adi(R,p4);
% calc for mobile twist (twist dagger ed)
ed=@(Adi,e)(Adi*e);
% about frame
edx=ed(Adix,twx);
edy=ed(Adiy,twy);
edz=ed(Adiz,twz);
edg=ed(Adig,twg);
edb=ed(Adib,twb);
eda=ed(Adia,twa);
O=zeros(6,1);
% about prop
  % prop1
  ed1x=ed(Adi1x,twx);
  ed1y=ed(Adi1y,twy);
  ed1z=ed(Adi1z,twz);
  ed1g=ed(Adi1g,twg);
  ed1b=ed(Adi1b,twb);
  ed1a=ed(Adi1a,twa);
  ed1=ed(Adi1,tw1);
  % prop 2
  ed2x=ed(Adi2x,twx);
  ed2y=ed(Adi2y,twy);
  ed2z=ed(Adi2z,twz);
  ed2g=ed(Adi2g,twg);
  ed2b=ed(Adi2b,twb);
  ed2a=ed(Adi2a,twa);
  ed2=ed(Adi2,tw2);
  %prop 3
  ed3x=ed(Adi3x,twx);
  ed3y=ed(Adi3y,twy);
  ed3z=ed(Adi3z,twz);
  ed3g=ed(Adi3g,twg);
  ed3b=ed(Adi3b,twb);
  ed3a=ed(Adi3a,twa);
  ed3=ed(Adi3,tw3);
  %prop 4
  ed4x=ed(Adi4x,twx);
  ed4y=ed(Adi4y,twy);
  ed4z=ed(Adi4z,twz);
  ed4g=ed(Adi4g,twg);
  ed4b=ed(Adi4b,twb);
  ed4a=ed(Adi4a,twa);
  ed4=ed(Adi4,tw4);
%% jacobian for frame
 Jo=[edx edy edz edg edb eda O zeros(6,3)];
%% jacobian for prop
    % prop 1
    J1=simplify([ed1x ed1y ed1z ed1g ed1b ed1a ed1 zeros(6,3)]);
    J2=simplify([ed2x ed2y ed2z ed2g ed2b ed2a ed2 zeros(6,3)]);
    J3=simplify([ed3x ed3y ed3z ed3g ed3b ed3a ed3 zeros(6,3)]);
    J4=simplify([ed4x ed4y ed4z ed4g ed4b ed4a ed4 zeros(6,3)]);
%% Moment of inertia matrix
% for frame
Moj=Jo'*Mo*Jo;
% for prop
Mij=@(J,M)(J.'*M*J);
M1j=Mij(J1,M1);
M2j=Mij(J2,M2);
M3j=Mij(J3,M3);
M4j=Mij(J4,M4);
% total inertia
M=simplify(4*Moj+M1j+M2j+M3j+M4j);
%% coriolis matrix C
% (a,b) element of C
C=sym(zeros(10,10));
D=[x,y,z,gamma,beta,alpha,th1,th2,th3,th4];
for i=1:10
    for j=1:10
        for k=1:10
            C(i,j)=C(i,j)+Cab(M,i,j,D,k);
        end
    end
end

%% function zone
% for Coriolis (the Tij)
function Tij=Cab(M,i,j,D,k)
    Tij=(diff(M(i,j),D(k))+diff(M(i,k),D(j))-diff(M(k,j),D(i)))/2;
end
