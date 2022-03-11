function [ZZZ,ZZ,ext,tuc]=Rectangle_biofilm_simulation
close all;
tic;
N=201; %Matrix width
Nrun=1000; %maybe 100 is enough
n=5; %neighbourhood size
nrep=25; %repression - not used
delta=0.1; %strength of selection
delta2=0.00001; %migration
repflag=0;
normflag=1; % normalization of fecundity by number of neighbors
plotflag=0;
expflag=0;

sr = strel('disk',1); %direct neighbours
sn = strel('disk',n); %social neighbourhood
snm=getnhood(sn);
srm=getnhood(sr);
srm(2,2)=0;
Tmax=500;

for kb=1:5
B=0.01*2^kb; %benefit parameter

for kc=1:5
C=B*0.2*kc; %cost parameter

for fec=1:5
% tic;

ext{kb,kc,fec}=zeros(1,Nrun);
ZZZ{kb,kc,fec}=zeros(Tmax,Nrun);

for mm=1:Nrun

A=zeros(5,N); %A is the matrix of genotypes
%A(1,:)=ceil(rand(1,N)*20);%
% A(1,:)=ceil(+(rand(1,N)>0.05))+1; %initial conditions
% A(2,:)=A(1,:);
A(1,:)=2;
A(1,N/2+0.5)=1;
% A(1,1:N)=1:N;
A(2,:)=A(1,:);

kk=1;
flag=1;
count=0;
Ac=[];

for i=1:Tmax
    % Fill holes in P - make internal viables, nonviable. A is still with holes
    Ao=A;
    P=A>0;
    Pp = padarray(P, [0 N], 'circular'); %duplicate the matrix periodically to both sides
    Pp = [ones(1,N*3);Pp];
    Ph=imfill(Pp,'holes');
    Ph=Ph(2:end,N+1:end-N);

    %matrix size control
    if max2(Ao(end-2:end,:))>0 % If last 3 rows are not empty
        Ao=[Ao;zeros(2,N)]; % Add two lines of zeros at the bottom
    end
    %Delete full rows (leave two full rows) (after filling holes)
    Im=min(Ph');
    Imn=find(Im==0,1,'first');
    if Imn>3
        Ac=[Ac;Ao(1:Imn-3,:)]; %%% 3 instead of 2 (to avoid overlap) %%% and change order with next row
        Ao=Ao(Imn-2:end,:);
    end

    A=Ao;   
    P=A>0;
    % Fill holes in P - make internal viables, nonviable. AGAIN (for size regulated matrix).
    % A is still with holes (and always stays this way)
    Pp = padarray(P, [0 N], 'circular'); %duplicate the matrix periodically to both sides
    Pp = [ones(1,N*3);Pp];
    Ph=imfill(Pp,'holes');
    Ph=Ph(2:end,N+1:end-N);

    Pn=imdilatepad(Ph,sr).*(1-Ph); %immediate unoccupied positions = new sites
    Pc=imdilatepad(Pn,sr).*Ph; %"living" occupied positions
    Pco=Pc;

    A1=(A==1 & Pc==1); % DOESN'T COUNT HOLES BECAUSE A=0 THERE, BUT THEY'RE NOT VIABLE ANYWAY...
    A2=(A==2 & Pc==1);
    
    Ic=find(Pc);
    In=find(Pn);
    
    G1=conv2pad(A1,snm); %finding how many active sites of genotype 1 around a focal active site
    G2=conv2pad(A2,snm); %finding how many active sites of genotype 2 around a focal active site
    G=G1./(G1+G2); %fraction of genotype #1 in nhood
    G(Pc==0)=0;
    
    %Repression % GO OVER 1
    if repflag==1
        Nlim=20;
        LimMat=(size(Ao,1):-1:-size(Ao,1))'*ones(1,nrep)/nrep;
        NutLim=(1-min(max(0,(Nlim/2+conv2pad(Pc,LimMat))),Nlim)/Nlim).*Pc;
    else
        NutLim=Pc*0+1;
    end
    
    if expflag==1
        % Expected relatedness   % GO OVER 2
        aa=ones(3,3);aa(2,2)=0;
        aN=conv2pad(Pc,aa);
        na=hist(aN(Ic),0:8);
    end
    
    % Omitting sites that will surely die - irrelevant because we filled the holes
    N1=(4-conv2pad(P,srm)); %number of empty neighbours (incorrect for top row but doesn't matter)
    Ne=N1.*Pc; %number of empty neighbours to active site
    Nn=(N1>0).*Pn; % relevant immidiate empty sites (won't "die" automatically because they have at least one empty neighbor)
    Nn0=(N1==0).*Pn; %immidiate empty sites that are holes (will die automatically upon birth)
    Nen0=imdilate(Nn0,srm).*Pc; %living sites that can grow into these holes

    %Fecundity functions
    switch fec
        case 1
            W = (1+B-C).*A1 + (1+B-C).*A2; % UNIFORM (ALL COOPERATORS / FACULTATORS)
        case 2
            W = (1+B*G-C).*A1 + (1+B*G).*A2; % COOPERATOR (#1) AND CHEATER (#2)
        case 3
            W = (1+B*(1-G)).*A1 + (1+B*(1-G)-C).*A2; % CHEATER (#1) AND COOPERATOR (#2)
        case 4
            W = (1+B-C*G).*A1 + (1+B-C*(1-G)).*A2; % LINEAR FACULTATIVE CHEATING
        case 5
            W = (1+B*(G.^2+(1-G).^2)-C*G).*A1 + (1+B*(G.^2+(1-G).^2)-C*(1-G)).*A2; % QS FACULTATIVE CHEATING
    end
    W=W.*NutLim;%*0+1;
    
    % For average calculations
    Birth=W; %before normalization - sums birth of active sites into all directions
    
    if normflag==1
    W(Ic)=W(Ic)./Ne(Ic); %normalize W to fecundity per neighbor empty site
    end
    Birth=Birth-W.*Nen0; %growth into a square surrounded on all four sides is effictivly like not growing
    
    %terminating
    I1=find(A(Ic)==1);
    ZZZ{kb,kc,fec}(i,mm)=length(I1)/length(Ic);
    if length(I1)==0
        flag=4;
        ext{kb,kc,fec}(mm)=i;
        break
    end
    
    %Printing
    if mod(i,50)==0 && plotflag==1
        count=count+1;
        figure(1)
        imrgb([Ac; A])
%       axis normal
%       maximize
        disp([mm i length(unique(A(Ic)))])
        %ZZZ(mm,count)=length(unique(A(Ic)));
        pause(0.05);

        I2=find(A(Ic)==2); %I2 is calculated only once every 50 generations
        
        f{mm}(count)=length(I1)./(length(I1)+length(I2));
        
        E=cov(G(Ic),2-A(Ic));
        Rt{mm}(count)=E(1,2)/E(2,2);
        %non-clonal relatedness
        Inc=find(G(Ic)~=1);
        Enc=cov(G(Ic(Inc)),2-A(Ic(Inc)));
        Rnc{mm}(count)=Enc(1,2)/Enc(2,2);
        %clonal fraction
        fc{mm}(count)=(length(Ic)-length(Inc))/length(find(A(Ic)==1)); %This is r_\infty
    end
        %more prints
%     if mod(i,500)==0
%         figure(10)
%         subplot(3,1,1)
%         plot(1:i-1,f{mm});%,1:i,f(1)+cumsum(dZp))
%         
%         subplot(3,1,2)
%         plot(1:i,Rt{mm},1:i,Rnc)
%         %flag=3;
%         subplot(3,1,3)
%         plot(1:i,fc)
%         pause(0.01)
% %        A(Ic)=ceil(rand(1,length(Ic))*20);
%     end
    
    %new sites
    I=find(Pn==1);
    CS=[1 0;-1 0;0 -1;0 1];
    Wntot=Pn*0;
    for j=1:4
        Wc{j}=Pn.*circshift(W,CS(j,:)); %measure 4 fitness values for a new site (for 4 possible births)
        % some of the births are impossible and are not counted
        g{j}=circshift(A,CS(j,:));
        Wntot=Wntot+Wc{j}; %sum of 4 fitness values = proportional to birth probability of a new site
    end
    
    Mat=[ones(length(I),1) delta*Wc{1}(I) delta*Wc{2}(I) delta*Wc{3}(I) delta*Wc{4}(I)];
    V=choose(Mat); % for any possible new site it finds the index from 1 to 5 that wins the rand lottery
    
%     PPP=A*0;
    for k=1:4
        A(I(find(V==k+1)))=g{k}(I(find(V==k+1))); %replace new site with the genotype that duplicated into it
%         PP=A*0;
%         PP(I(find(V==k+1)))=1;
%         PP=PP.*Nn;
%         PPP=PPP+circshift(PP,-CS(k,:));
    end
    
    % MIGRATION
    CS1=[1 0;-1 0;0 -1;0 1;1 1;1 -1;-1 1;-1 -1];
    DiffFlag=0;
    switch DiffFlag
        case 1
            Pca = padarray(Pc, [0 1], 'circular');
            Aa = padarray(A, [0 1], 'circular');
            Ica=find(Pca);
            DiffIndex=(rand(1,length(Ica))<delta2).*ceil(rand(1,length(Ica))*8);
            [Icc,Jcc]=find(Pca);
            for k=1:8
                Id=find(DiffIndex==k);
                Id=Id(find(Jcc(Id)>1 & Jcc(Id)<N+2));
                Idd = sub2ind(size(Pca), Icc(Id)+CS1(k,1), Jcc(Id)+CS1(k,2));
                Ad=Aa(Idd).*Pca(Idd);
                Idn=find(Ad~=0);
                VV=Aa(Idd(Idn));
                Aa(Idd(Idn))=Aa(Ica(Id(Idn)));
                Aa(Ica(Id(Idn)))=VV;
            end
            A=Aa(:,2:N+1);
        case 2 %global
            DiffIndex=find(rand(1,length(Ic))<delta2);
            Ndiff=length(DiffIndex);
            A(Ic(DiffIndex(randperm(Ndiff))))=A(Ic(DiffIndex));
    end

    %expected changes
%     Death=conv2pad(Wntot,srm).*(Ne==1);
%     Wt=1+delta*(Birth(Ic)-Death(Ic));
%     dZp(i)=mean(((2-A(Ic))-mean(2-A(Ic))).*(Wt-mean(Wt)))/mean(Wt);
%     
%     P=A>0;
%     Pp = padarray(P, [0 3], 'circular');
%     Pp= [ones(1,N+6);Pp];
%     Ph=imfill(Pp,'holes'); %fill holes in P - make internal viables, nonviable.
%     Ph=Ph(2:end,4:end-3);
    
%     Pco=Pc;
%     Pn=imdilatepad(Ph,sr).*(1-Ph);
%     Pc=imdilatepad(Pn,sr).*Ph;
    
%     ActualFittness=Pco+PPP-(Pco==1 & Pc==0);
%     dft(i)=(mean(ActualFittness(Ic).*(2-A(Ic)))-mean(ActualFittness(Ic)).*mean(2-A(Ic)))/mean(ActualFittness(Ic));
%     A1=(A==1 & Pc==1);
%     A2=(A==2 & Pc==1);
%     f(i)=sum(sum(A1))/(sum(sum(A2))+sum(sum(A1)));
%     
%     nI(i)=length(Ic);

%     if i>1 & abs(dft(i)-(f(i)-f(i-1)))>1e-10
%     end

%     if flag==3
%         tr(kk)=f(i);
%         A=Ao;
%         kk=kk+1;
%     end
end

if flag==4
    ZZZ{kb,kc,fec}(i+1:Tmax,mm)=0;
    %ZZZ(mm)=length(I1)/length(Ic);
    flag=0;
end

end
ZZ{kb,kc,fec}=mean(ZZZ{kb,kc,fec},2);
length(ext{kb,kc,fec}==0);
tuc{kb,kc,fec}=toc;
end

% figure;
% for fec=1:5
%     subplot(2,1,1); plot(ZZ{kb,kc,fec}); hold on;
%     ylim([0 1]);
%     subplot(2,1,2); histogram(ext{kb,kc,fec},Tmax); hold on;
%     xlim([0 Tmax]); ylim([0 Nrun]);
% end
end
save(strcat('run',num2str(kb),'.mat'));
end
return

function dP=imdilatepad(P,s)
sn=getnhood(s);
snz=(size(sn)+1)/2; % actually this is just [2,2]
Pz=size(P);
Pp = padarray(P, [0 snz(2)], 'circular'); % add 2 columns before and after the vacancy matrix with periodical boundry cond.
dPp=imdilate(Pp,s); % add a "row" of ones at the interface (occupy a "row" of empty places)
dP=dPp((1:Pz(1)),snz(2)+(1:Pz(2))); % return to the original size matrix (without the added boundries)

function dP=conv2pad(P,s)
if ~isnumeric(+s)
    sn=getnhood(s);
else
    sn=s;
end
snz=(size(sn)+1)/2;
Pz=size(P);
Pp = padarray(P, [0 snz(1)], 'circular');
dP=conv2(+Pp,+sn,'same');
dP=dP(:,snz(1)+1:end-snz(1));

function O=choose(Mat)
M=Mat';
[m,n]=size(M);
Mc=sign(cumsum(M)-ones(m,1)*(rand(1,n).*sum(M)));
Mt=[-ones(1,n);Mc];
[O,~]=find(diff(Mt)); % finds the index from 1 to 5 that wins the rand lottery

% function varargout=meanPlot(x,y,x0)