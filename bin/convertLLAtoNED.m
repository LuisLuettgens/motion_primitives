function [x,y,z]=convertLLAtoNED(long,lat,h,long0,lat0,h0)
    
    %long0=53.2;
    %lat0=8.5;
    %long0=12.10074142;
    %lat0=54.17057475;
    %long0=12.10981855183098;
    %lat0=54.16653821183098;
    %h0=0;

    R_ea=6378137.0;%[m]
    f=1/298.257223563;
    R_eb=R_ea*(1-f);%[m]
    e=sqrt(R_ea^2-R_eb^2)/R_ea;
    
    N_e=R_ea/sqrt(1-e^2*sind(lat)^2);
    N_e0=R_ea/sqrt(1-e^2*sind(lat0)^2);
    
    P_e=[(N_e+h)*cosd(lat)*cosd(long); (N_e+h)*cosd(lat)*sind(long); (N_e*(1-e^2)+h)*sind(lat)];
    P_e0=[(N_e0+h0)*cosd(lat0)*cosd(long0); (N_e0+h0)*cosd(lat0)*sind(long0); (N_e0*(1-e^2)+h0)*sind(lat0)];
    
    R_ne=[-sind(lat0)*cosd(long0) -sind(lat0)*sind(long0) cosd(lat0); -sind(long0) cosd(long0) 0.; -cosd(lat0)*cosd(long0) -cosd(lat0)*sind(long0) -sind(lat0)];
    
    P_n=R_ne*(P_e-P_e0);
    
    x=P_n(1);
    y=P_n(2);
    z=P_n(3);

end