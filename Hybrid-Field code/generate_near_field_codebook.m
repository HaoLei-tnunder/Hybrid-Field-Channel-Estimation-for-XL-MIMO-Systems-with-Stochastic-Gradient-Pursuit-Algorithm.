function [codebook,record] = generate_near_field_codebook(N1,N2,d,P1,P2,Delta)

N=N1*N2;

Xmax1=P1(1); Xmin1=P1(2); Ymax1=P1(3); Ymin1=P1(4); Zmax1=P1(5); Zmin1=P1(6);
Xmax2=P2(1); Xmin2=P2(2); Ymax2=P2(3); Ymin2=P2(4); Zmax2=P2(5); Zmin2=P2(6);
Xdelta1=Delta(1);Ydelta1=Delta(2);Zdelta1=Delta(3);Xdelta2=Delta(4);Ydelta2=Delta(5);Zdelta2=Delta(6);

% Xgrid1 = linspace(Xmin1,Xmax1,Xnum1); Ygrid1 = linspace(Ymin1,Ymax1,Ynum1); Zgrid1 = linspace(Zmin1,Zmax1,Znum1);
% Xgrid2 = linspace(Xmin2,Xmax2,Xnum2); Ygrid2 = linspace(Ymin2,Ymax2,Ynum2); Zgrid2 = linspace(Zmin2,Zmax2,Znum2);

Xgrid1=[Xmin1:Xdelta1:Xmax1]; Ygrid1=[Ymin1:Ydelta1:Ymax1]; Zgrid1=[Zmin1:Zdelta1:Zmax1]; 
Xgrid2=[Xmin2:Xdelta2:Xmax2]; Ygrid2=[Ymin2:Ydelta2:Ymax2]; Zgrid2=[Zmin2:Zdelta2:Zmax2]; 

Xnum1=length(Xgrid1); Ynum1=length(Ygrid1); Znum1=length(Zgrid1);
Xnum2=length(Xgrid2); Ynum2=length(Ygrid2); Znum2=length(Zgrid2);

record=zeros(Xnum1*Ynum1*Znum1*Xnum2*Ynum2*Znum2,6);
codebook = zeros(Xnum1*Ynum1*Znum1*Xnum2*Ynum2*Znum2,N);
i=1;
for x1=Xgrid1
    for y1=Ygrid1
        for z1=Zgrid1
            for x2=Xgrid2
                for y2=Ygrid2
                    for z2=Zgrid2
                        a = zeros(1,N);
                        for n1=1:N1
                            for n2=1:N2
                                a((n1-1)*N2+n2) = exp(1j*2*pi*(sqrt((x1-(n1-1-(N1-1)/2)*d)^2+(z1-(n2-1-(N2-1)/2)*d)^2+y1^2)+sqrt((x2-(n1-1-(N1-1)/2)*d)^2+(z2-(n2-1-(N2-1)/2)*d)^2+y2^2)));
                            end
                        end
%                         a=a/sqrt(N);
                        codebook(i,:)=a;
                        record(i,:)=[x1,y1,z1,x2,y2,z2];
                        i=i+1;
                    end
                end
            end
        end
    end
end

[codebook,index]=unique(codebook,'row');
record=record(index,:);
end

