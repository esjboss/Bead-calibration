 
function [VSF]=VSF_beads(D0, dD, lambda0, dlambda, angle)

% function computing the phase function of polystiren calibration beads
% assuming both them and the wavelength-range are both normally distributed
% D0-mean size [micron]
% dD - standard deviation of size [micron]
% lambda - wavelength in air [nm]
% dlambda - standard deviation of wavelength distribution[nm]
% angle - output angles [degrees] 
% if dD == 0, only nominal size is used.
% if dlambda == 0, only nominal wavelength is used.

%%% used in all loops
nang=100*3+1; %number of angles for mie computation
dang=pi/2/(nang-1); %angular resolution
nan=2*nang-1; %number of angles for computation 
ang=[0:dang:pi]*180/pi;

lambda0=lambda0/1000;%convert to microns
dlambda=dlambda/1000;

if dD >0  %size weighing function parameters
    min_D=D0-3*dD;
    max_D=D0+3*dD;
    N=100;
    D_inc=(max_D-min_D)/N;
    d=[min_D:D_inc:max_D];
    comp_d=N+1; %number of wavelengths computations
    D_func=exp(-0.5*((d-D0)/dD).^2);
end

if dlambda >0 %wavelength weighing function parameters
    min_wl=lambda0-3*dlambda;
    max_wl=lambda0+3*dlambda;
    N=100;
    wl_inc=(max_wl-min_wl)/N;
    wl=[min_wl:wl_inc:max_wl];
    comp_wl=N+1; %number of wavelengths computations
    W_func=exp(-0.5*((wl-lambda0)/dlambda).^2);
end

%input for index of refraction of water
S=0;
T=20;

if (dD == 0  & dlambda == 0)
    lambda0
    [n, nm]=IoR(lambda0,S,T);
    rho=pi*D0/lambda0*nm;
    [S1 S2 Qb Qc Qback]=bhmie(rho,n,nang);
    S11=0.5*((abs(S1)).^2 + (abs(S2)).^2);
    S11_rad=S11.*sin(ang'/180*pi)*2*pi*dang;
    S11_int=3/8*(S11_rad(1)+S11_rad(nan))+7/6*(S11_rad(2)+S11_rad(nan-1))+23/24*(S11_rad(3)+S11_rad(nan-2));
    S11_int=(S11_int+sum(S11_rad(4:nan-3)));
    S11=S11/S11_int; %phase function
    S11_ang=interp1(ang,S11,angle,'linear');
    VSF=S11_ang*Qb*pi*D0^2/4;
    
elseif (dD == 0 & dlambda>0)
    for j=1:comp_wl
        lambda(j)=min_wl+(j-1)*wl_inc;
        weight_wl(j)=interp1(wl,W_func,lambda(j),'pchip');
        [n, nm]=IoR(lambda(j),S,T);
        rho=pi*D0/lambda(j)*nm; %rho
        [S1 S2 Qb Qc Qback]=bhmie(rho,n,nang);
        S11=0.5*((abs(S1)).^2 + (abs(S2)).^2);
        S11_rad=S11.*sin(ang'/180*pi)*2*pi*dang;
        S11_int=3/8*(S11_rad(1)+S11_rad(nan))+7/6*(S11_rad(2)+S11_rad(nan-1))+23/24*(S11_rad(3)+S11_rad(nan-2));
        S11_int=(S11_int+sum(S11_rad(4:nan-3)));
        S11=S11/S11_int; %phase function
        factor=Qb*pi*D0^2/4;
        B(j,:)=weight_wl(j)*S11*factor*wl_inc;
    end
    %wavelength integration
    int_wl=3/8*(weight_wl(1)+weight_wl(comp_wl))+7/6*(weight_wl(2)+weight_wl(comp_wl-1))+23/24*(weight_wl(3)+weight_wl(comp_wl-2));
    int_wl=(int_wl+sum(weight_wl(4:comp_wl-3)))*wl_inc;
    
    B_int=3/8*(B(1,:)+B(comp_wl,:))+7/6*(B(2,:)+B(comp_wl-1,:))+23/24*(B(3,:)+B(comp_wl-2,:));
    VSF_=(B_int+sum(B(4:comp_wl-3,:)))/int_wl;
    VSF=interp1(ang,VSF_,angle,'linear');
    
    elseif (dD > 0 & dlambda == 0)
    [n, nm]=IoR(lambda0,S,T);
    for j=1:comp_d
        D(j)=min_D+(j-1)*D_inc;
        weight_D(j)=interp1(d,D_func,D(j),'pchip');        
        rho=pi*D(j)/lambda0*nm; %rho
        [S1 S2 Qb Qc Qback]=bhmie(rho,n,nang);
        S11=0.5*((abs(S1)).^2 + (abs(S2)).^2);
        S11_rad=S11.*sin(ang'/180*pi)*2*pi*dang;
        S11_int=3/8*(S11_rad(1)+S11_rad(nan))+7/6*(S11_rad(2)+S11_rad(nan-1))+23/24*(S11_rad(3)+S11_rad(nan-2));
        S11_int=(S11_int+sum(S11_rad(4:nan-3)));
        S11=S11/S11_int; %phase function
        factor=Qb*pi*D(j)^2/4;
        B(j,:)=weight_D(j)*S11*factor*D_inc;
    end
    %wavelength integration
    int_D=3/8*(weight_D(1)+weight_D(comp_d))+7/6*(weight_D(2)+weight_D(comp_d-1))+23/24*(weight_D(3)+weight_D(comp_d-2));
    int_D=(int_D+sum(weight_D(4:comp_d-3)))*D_inc;
    
    B_int=3/8*(B(1,:)+B(comp_d,:))+7/6*(B(2,:)+B(comp_d-1,:))+23/24*(B(3,:)+B(comp_d-2,:));
    VSF_=(B_int+sum(B(4:comp_d-3,:)))/int_D;
    VSF=interp1(ang,VSF_,angle,'linear');


%full monty
elseif (dlambda > 0 & dD>0)
    for j=1:comp_d
        lambda(j)=min_wl+(j-1)*wl_inc;
        weight_wl(j)=interp1(wl,W_func,lambda(j),'pchip');
        for m=1:comp_d
            D(m)=min_D+(m-1)*D_inc;
            weight_D(m)=interp1(d,D_func,D(m),'pchip');        
            [n, nm]=IoR(lambda(j),S,T);
            rho=pi*D(m)/lambda(j)*nm; %rho
            [S1 S2 Qb Qc Qback]=bhmie(rho,n,nang);
            S11=0.5*((abs(S1)).^2 + (abs(S2)).^2);
            S11_rad=S11.*sin(ang'/180*pi)*2*pi*dang;
            S11_int=3/8*(S11_rad(1)+S11_rad(nan))+7/6*(S11_rad(2)+S11_rad(nan-1))+23/24*(S11_rad(3)+S11_rad(nan-2));
            S11_int=(S11_int+sum(S11_rad(4:nan-3)));
            S11=S11/S11_int; %phase function
            factor=Qb*pi*D(j)^2/4;
            
            B(m,:)=factor*weight_D(m)*weight_wl(j)*S11*D_inc*wl_inc;
        end
        %size integration
        B_int1(j,:)=3/8*(B(1,:)+B(comp_d,:))+7/6*(B(2,:)+B(comp_d-1,:))+23/24*(B(3,:)+B(comp_d-2,:));
        B_int1(j,:)=B_int1(j,:)+sum(B(4:comp_d-3,:));
    end
    %wavelength integration
    int_wl=3/8*(weight_wl(1)+weight_wl(comp_wl))+7/6*(weight_wl(2)+weight_wl(comp_wl-1))+23/24*(weight_wl(3)+weight_wl(comp_wl-2));
    int_wl=(int_wl+sum(weight_wl(4:comp_wl-3)))*wl_inc;
    
    int_D=3/8*(weight_D(1)+weight_D(comp_d))+7/6*(weight_D(2)+weight_D(comp_d-1))+23/24*(weight_D(3)+weight_D(comp_d-2));
    int_D=(int_D+sum(weight_D(4:comp_d-3)))*D_inc;

    
    B_int2=3/8*(B_int1(1,:)+B_int1(comp_wl,:))+7/6*(B_int1(2,:)+B_int1(comp_wl-1,:))+23/24*(B_int1(3,:)+B_int1(comp_wl-2,:));
    B_int2=(B_int2+sum(B_int1(4:comp_wl-3,:)))/int_wl/int_D;
    VSF=interp1(ang,B_int2,angle,'linear');
end
return

function [n, nm]=IoR(wl,S,T)

np=1.5718 + 0.008412/(wl^2) + 0.000235/(wl^4); %Jones et al.
ni=0.0003;  %imaginary part of index of refraction
n0=1.31405; n1=1.779e-4; n2=-1.05e-6; n3=1.6e-8; n4=-2.02e-6; n5=15.868; n6=0.01155; n7=-0.00423; n8=-4382; n9=1.1455e6; %Quan and Fry, 1995
nm=n0+(n1+n2*T+n3*T^2)*S+n4*T^2+(n5+n6*S+n7*T)/(wl*1000)+n8/(wl*1000)^2+n9/(wl*1000)^3; %checked
n=(np+ni*sqrt(-1))/nm;

return


