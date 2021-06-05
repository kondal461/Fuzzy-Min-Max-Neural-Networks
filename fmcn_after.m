%% clear command window

clc;
clear all;
close all;


%% Input



load Iris_dataset
p=Iris_dataset;
x=p(:,1:4);
x_label=p(:,end);         %1*irisTargets(1,:)+2*irisTargets(2,:)+3*irisTargets(3,:);
n_class=max(x_label);
%x_label=x_label';
%x= y./max(max(y));

[Train, Test] = crossvalind('HoldOut', x_label, 0.5);

TrainingSample=x(Train,:);
TrainingLabel=x_label(Train,:);
TestingSample=x(Test,:);
TestingLabel=x_label(Test,:);

input1=TrainingSample((TrainingLabel==1),:);
input2=TrainingSample((TrainingLabel==2),:);
input3=TrainingSample((TrainingLabel==3),:);

test1=TestingSample((TestingLabel==1),:);
test2=TestingSample((TestingLabel==2),:);
test3=TestingSample((TestingLabel==3),:);


nFeature=size(x,2);
nTrainData=size(input1,1);


%% hyperbox Creation & Expansion


theta=0.15
gamma=4;

l_theta=length(theta);
tic
for iTrail=1:200
for th=1:l_theta
    
    
    Vj=zeros(nTrainData,nFeature);
    Wj=zeros(nTrainData,nFeature);
    Vj_new=zeros(nTrainData,nFeature);
    Wj_new=zeros(nTrainData,nFeature);
    x=dataset();
    Ix=1;
    
    Vk=zeros(nTrainData,nFeature);
    Wk=zeros(nTrainData,nFeature);
    Vk_new=zeros(nTrainData,nFeature);
    Wk_new=zeros(nTrainData,nFeature);
    y=dataset();
    Iy=1;
    
    Vl=zeros(nTrainData,nFeature);
    Wl=zeros(nTrainData,nFeature);
    Vl_new=zeros(nTrainData,nFeature);
    Wl_new=zeros(nTrainData,nFeature);
    z=dataset();
    Iz=1;
    
    for j=1
        Vj(j,:)=input1(j,:);
        Wj(j,:)=input1(j,:);
        Vk(j,:)=input2(j,:);
        Wk(j,:)=input2(j,:);
        Vl(j,:)=input3(j,:);
        Wl(j,:)=input3(j,:);
        
        x.hyperbox{j}=[Vj(j,:);Wj(j,:)];
        y.hyperbox{j}=[Vk(j,:);Wk(j,:)];
        z.hyperbox{j}=[Vl(j,:);Wl(j,:)];
    end
    for j=2:size(input1,1)
        x1_hyp=x.hyperbox;
        y1_hyp=y.hyperbox;
        z1_hyp=z.hyperbox;
        
        x1=[x1_hyp{Ix}(1,:);input1(j,:)];
        x2=[x1_hyp{Ix}(2,:);input1(j,:)];
        p=HB_expansion(x2,x1);
        
        x3=[y1_hyp{Iy}(1,:);input2(j,:)];
        x4=[y1_hyp{Iy}(2,:);input2(j,:)];
        p1=HB_expansion(x4,x3);
        
        x5=[z1_hyp{Iz}(1,:);input3(j,:)];
        x6=[z1_hyp{Iz}(2,:);input3(j,:)];
        p2=HB_expansion(x6,x5);
        
        if(nFeature*theta(th)>=p)&&(isolation1(x1_hyp(Ix),y1_hyp(1:Iy),nFeature))&&(isolation1(x1_hyp(Ix),z1_hyp(1:Iz),nFeature))
            Vj_new(j,:)=min(x1,[],1);
            Wj_new(j,:)=max(x2,[],1);
            Vj(j,:)=Vj_new(j,:);
            Wj(j,:)=Wj_new(j,:);
            if(j == nTrainData)
                x.hyperbox{Ix}=[Vj_new(j,:);Wj_new(j,:)];
                x1_hyp=x.hyperbox;
            end
            x.hyperbox{Ix}=[Vj(j,:);Wj(j,:)];
            x1_hyp=x.hyperbox;
        else
            Vj(j,:)=input1(j,:);
            Wj(j,:)=input1(j,:);
            Ix=Ix+1;
            x.hyperbox{Ix}=[Vj(j,:);Wj(j,:)];
            x1_hyp=x.hyperbox;
        end
        
        if(nFeature*theta(th)>=p1)&&(isolation1(y1_hyp(Iy),x1_hyp(1:Ix),nFeature))&&(isolation1(y1_hyp(Iy),z1_hyp(1:Iz),nFeature))
            Vk_new(j,:)=min(x3,[],1);
            Wk_new(j,:)=max(x4,[],1);
            Vk(j,:)=Vk_new(j,:);
            Wk(j,:)=Wk_new(j,:);
            if(j==nTrainData)
                y.hyperbox{Iy}=[Vk_new(j,:);Wk_new(j,:)];
                y1_hyp=y.hyperbox;
            end
            y.hyperbox{Iy}=[Vk(j,:);Wk(j,:)];
            y1_hyp=y.hyperbox;
        else
            Vk(j,:)=input2(j,:);
            Wk(j,:)=input2(j,:);
            Iy=Iy+1;
            y.hyperbox{Iy}=[Vk(j,:);Wk(j,:)];
            y1_hyp=y.hyperbox;
        end
        
        if(nFeature*theta(th)>=p2)&&(isolation1(z1_hyp(Iz),x1_hyp(1:Ix),nFeature))&&(isolation1(z1_hyp(Iz),y1_hyp(1:Iy),nFeature))
            Vl_new(j,:)=min(x5,[],1);
            Wl_new(j,:)=max(x6,[],1);
            Vl(j,:)=Vl_new(j,:);
            Wl(j,:)=Wl_new(j,:);
            if(j==nTrainData)
                z.hyperbox{Iz}=[Vl_new(j,:);Wl_new(j,:)];
                z1_hyp=z.hyperbox;
            end
            z.hyperbox{Iz}=[Vl(j,:);Wl(j,:)];
            z1_hyp=z.hyperbox;
        else
            Vl(j,:)=input3(j,:);
            Wl(j,:)=input3(j,:);
            Iz=Iz+1;
            z.hyperbox{Iz}=[Vl(j,:);Wl(j,:)];
            z1_hyp=z.hyperbox;
        end
        
    end
    x1=dataset();
    y1=dataset();
    z1=dataset();
    
    for i=1:size(x1_hyp,2)
        if(x1_hyp{i}~=inf(2,nFeature))
            x1.hyperbox{i}=x1_hyp{i};
        end
    end
    for i=1:size(y1_hyp,2)
        if(y1_hyp{i}~=inf(2,nFeature))
            y1.hyperbox{i}=y1_hyp{i};
        end
    end
    for i=1:size(z1_hyp,2)
        if(z1_hyp{i}~=inf(2,nFeature))
            z1.hyperbox{i}=z1_hyp{i};
        end
    end

    
    x_hyp=x1.hyperbox;  % No of hyperbox created in class-1
    y_hyp=y1.hyperbox;  % No of hyperbox created in class-2
    z_hyp=z1.hyperbox;  % No of hyperbox created in class-3
    
%     fprintf('x_hyp\n');
%     disp(x_hyp);
%     fprintf('y_hyp\n');
%     disp(y_hyp);
%     fprintf('z_hyp\n');
%     disp(z_hyp);
%     
    
    
    
    
    %% containment & overlap compensate neuron
   
    Xc=dataset();
    Xo=dataset();
    a=1;
    
    a1=1;
    for i=1:size(x_hyp,2)
        for j=1:size(y_hyp,2)
            if(isolation_test(x_hyp(i),y_hyp(j)))
            elseif(Cont_Comp_Neuron(x_hyp(i),y_hyp(j)))
                Vc=max(y_hyp{j}(1,:),x_hyp{i}(1,:));
                Wc=min(y_hyp{j}(2,:),x_hyp{i}(2,:));
                Xc.containment{a}=[Vc;Wc];
                X_Cont=Xc.containment;
                a=a+1;
            else
                if(x_hyp{i}(1,:)~=x_hyp{i}(2,:))&(y_hyp{j}(1,:)~=y_hyp{j}(2,:))
                    Vo=max(y_hyp{j}(1,:),x_hyp{i}(1,:));
                    Wo=min(y_hyp{j}(2,:),x_hyp{i}(2,:));
                    Xo.overlap{a1}=[Vo;Wo];
                    X_Over=Xo.overlap;
                    a1=a1+1;
                end
            end
        end
    end
    
    for i=1:size(x_hyp,2)
        for j=1:size(z_hyp,2)
            if(isolation_test(x_hyp(i),z_hyp(j)))
            elseif(Cont_Comp_Neuron(x_hyp(i),z_hyp(j)))
                Vc=max(z_hyp{j}(1,:),x_hyp{i}(1,:));
                Wc=min(z_hyp{j}(2,:),x_hyp{i}(2,:));
                Xc.containment{a}=[Vc;Wc];
                X_Cont=Xc.containment;
                a=a+1;
            else
                if(x_hyp{i}(1,:)~=x_hyp{i}(2,:))&(z_hyp{j}(1,:)~=z_hyp{j}(2,:))
                    Vo=max(z_hyp{j}(1,:),x_hyp{i}(1,:));
                    Wo=min(z_hyp{j}(2,:),x_hyp{i}(2,:));
                    Xo.overlap{a1}=[Vo;Wo];
                    X_Over=Xo.overlap;
                    a1=a1+1;
                end
            end
        end
    end
    
    for i=1:size(y_hyp,2)
        for j=1:size(z_hyp,2)
            if(isolation_test(y_hyp(i),z_hyp(j)))
            elseif(Cont_Comp_Neuron(y_hyp(i),z_hyp(j)))
                Vc=max(z_hyp{j}(1,:),y_hyp{i}(1,:));
                Wc=min(z_hyp{j}(2,:),y_hyp{i}(2,:));
                Xc.containment{a}=[Vc;Wc];
                X_Cont=Xc.containment;
                a=a+1;
            else
                if(y_hyp{i}(1,:)~=y_hyp{i}(2,:))&(z_hyp{j}(1,:)~=z_hyp{j}(2,:))
                    Vo=max(z_hyp{j}(1,:),y_hyp{i}(1,:));
                    Wo=min(z_hyp{j}(2,:),y_hyp{i}(2,:));
                    Xo.overlap{a1}=[Vo;Wo];
                    X_Over=Xo.overlap;
                    a1=a1+1;
                end
            end
        end
    end
    
    
    tf = isempty(Xc);
    if(tf==1)
        X_Cont{1}=[zeros(1,nFeature);zeros(1,nFeature)];
    end
    
    tf1 = isempty(Xo);
    if(tf1==1)
        X_Over{1}=[zeros(1,nFeature);zeros(1,nFeature)];
    end
    
    X_Cont; % Containment Compensation Neuron
    X_Over; % Overlap Compensation Neuron
    
    %% Learning or Training
    % DESCRIPTIVE TEXT
    test=[test1;test2;test3];
    n_test=size(test,1);
    n_test1=size(test1,1);
    n_test2=size(test2,1);
    n_test3=size(test3,1);
    count1=zeros(1,n_test);
    count2=zeros(1,n_test);
    count3=zeros(1,n_test);
    count4=zeros(1,n_test);
    i1=1;
    for r=1:size(test,1)
        
        Ah=test(r,:);
        
        
        %membership func for class-1
        
        B_c1=zeros(1,size(x_hyp,2));
        for j=1:size(x_hyp,2)
            B_c1(j)=MembershipFunc(Ah,x_hyp{j}(1,:),x_hyp{j}(2,:),gamma);
        end
        
        %membership func for class-2
        
        B_c2=zeros(1,size(y_hyp,2));
        for j=1:size(y_hyp,2)
            B_c2(j)=MembershipFunc(Ah,y_hyp{j}(1,:),y_hyp{j}(2,:),gamma);
        end
        
        %membership func for class-3
        B_c3=zeros(1,size(z_hyp,2));
        for j=1:size(z_hyp,2)
            B_c3(j)=MembershipFunc(Ah,z_hyp{j}(1,:),z_hyp{j}(2,:),gamma);
        end
        
        
        bj=B_c1;
        bk=B_c2;
        bl=B_c3;
        
        disp(['B1=',num2str(B_c1)]);
        disp(['B2=',num2str(B_c2)]);
        disp(['B3=',num2str(B_c3)]);
        
        
        
        %% Activation func of containment compensation
        % DESCRIPTIVE TEXT
        
        Xc1=dataset();
        Xc2=dataset();
        Xc3=dataset();
        
        e_ccn_c1=0;
        e_ccn_c2=0;
        e_ccn_c3=0;
        
        for k=1:size(X_Cont,2)
            if ( Ah(1:end)<=X_Cont{k}(2,1:end) & Ah(1:end)>=X_Cont{k}(1,1:end))
                
                m1=1;
                m2=1;
                m3=1;
                Xc1.contain1{m1}=[zeros(1,nFeature);zeros(1,nFeature)];
                for i=1:size(bj,2)           %select hyperbox those are containment
                    if (bj(i)==1)
                        g1=x_hyp{i};
                        Xc1.contain1{m1}=g1;
                        m1=m1+1;
                    end
                end
                
                Xc2.contain2{m2}=[zeros(1,nFeature);zeros(1,nFeature)];
                for i=1:size(bk,2)             %select hyperbox those are containment
                    if (bk(i)==1)
                        g2=y_hyp{i};
                        Xc2.contain2{m2}=g2;
                        m2=m2+1;
                    end
                end
                
                Xc3.contain3{m3}=[zeros(1,nFeature);zeros(1,nFeature)];
                for i=1:size(bl,2)             %select hyperbox those are containment
                    if (bl(i)==1)
                        g3=z_hyp{i};
                        Xc3.contain3{m3}=g3;
                        m3=m3+1;
                    end
                end
                
                Xc_1=Xc1.contain1;
                Xc_2=Xc2.contain2;
                Xc_3=Xc3.contain3;
                
                for i=1:size(Xc_1,2)     %the hyperbox where the data falls
                    for j=1:size(Xc_2,2)
                        if (max(Xc_1{i}(1,:),Xc_2{j}(1,:))==X_Cont{k}(1,1:end))&(min(Xc_1{i}(2,:),Xc_2{j}(2,:))==X_Cont{k}(2,1:end))
                            g4=Xc_1{i};
                            g5=Xc_2{j};
                            [ej_ccn]=CCN_ActivationFun(Ah,X_Cont{k}(1,:),X_Cont{k}(2,:),gamma);
                            if(containment_check1(g4,g5))
                                e_ccn_c1=ej_ccn;
                                e_ccn_c2=0;
                            elseif(containment_check2(g4,g5))
                                e_ccn_c1=ej_ccn;
                                e_ccn_c2=ej_ccn;
                            else
                                e_ccn_c1=0;
                                e_ccn_c2=ej_ccn;
                                
                            end
                        end
                    end
                end
                
                for i=1:size(Xc_1,2)     %the hyperbox where the data falls
                    for j=1:size(Xc_3,2)
                        if (max(Xc_1{i}(1,:),Xc_3{j}(1,:))==X_Cont{k}(1,1:end))&(min(Xc_1{i}(2,:),Xc_3{j}(2,:))==X_Cont{k}(2,1:end))
                            g4=Xc_1{i};
                            g5=Xc_3{j};
                            [ej_ccn]=CCN_ActivationFun(Ah,X_Cont{k}(1,:),X_Cont{k}(2,:),gamma);
                            if(containment_check1(g4,g5))
                                e_ccn_c1=ej_ccn;
                                e_ccn_c3=0;
                            elseif(containment_check2(g4,g5))
                                e_ccn_c1=ej_ccn;
                                e_ccn_c3=ej_ccn;
                            else
                                e_ccn_c1=0;
                                e_ccn_c3=ej_ccn;
                                
                            end
                        end
                    end
                end
                for i=1:size(Xc_2,2)     %the hyperbox where the data falls
                    for j=1:size(Xc_3,2)
                        if (max(Xc_2{i}(1,:),Xc_3{j}(1,:))==X_Cont{k}(1,1:end))&(min(Xc_2{i}(2,:),Xc_3{j}(2,:))==X_Cont{k}(2,1:end))
                            g4=Xc_2{i};
                            g5=Xc_3{j};
                            [ej_ccn]=CCN_ActivationFun(Ah,X_Cont{k}(1,:),X_Cont{k}(2,:),gamma);
                            if(containment_check1(g4,g5))
                                e_ccn_c2=ej_ccn;
                                e_ccn_c3=0;
                            elseif(containment_check2(g4,g5))
                                e_ccn_c2=ej_ccn;
                                e_ccn_c3=ej_ccn;
                            else
                                e_ccn_c2=0;
                                e_ccn_c3=ej_ccn;
                                
                            end
                        end
                    end
                end
                
            end
        end
        disp(['e_ccn_c1=',num2str(e_ccn_c1)]);
        disp(['e_ccn_c2=',num2str(e_ccn_c2)]);
        disp(['e_ccn_c3=',num2str(e_ccn_c3)]);
        
        
        
        
        %% Activation function of OCN
       
        Xo1=dataset();
        Xo2=dataset();
        Xo3=dataset();
        
        dp1_c1=0;
        dp2_c2=0;
        dp3_c3=0;
        for k=1:size(X_Over,2)
            if ( Ah(1:end)<=X_Over{k}(2,1:end) & Ah(1:end)>=X_Over{k}(1,1:end))
                m1=1;
                m2=1;
                m3=1;
                Xo1.over1{m1}=[zeros(1,n_fea);zeros(1,nFeature)];
                for i=1:size(bj,2)
                    if (bj(i)==1)
                        g6=x_hyp{i};
                        Xo1.over1{m1}=g6;
                        m1=m1+1;
                    end
                end
                
                Xo2.over2{m2}=[zeros(1,n_fea);zeros(1,n_fea)];
                for i=1:size(bk,2)
                    if (bk(i)==1)
                        g7=y_hyp{i};
                        Xo2.over2{m2}=g7;
                        m2=m2+1;
                    end
                end
                
                Xo3.over3{m3}=[zeros(1,n_fea);zeros(1,n_fea)];
                for i=1:size(bl,2)
                    if (bl(i)==1)
                        g8=z_hyp{i};
                        Xo3.over3{m3}=g8;
                        m3=m3+1;
                    end
                end
                
                Xo_1=Xo1.over1;
                Xo_2=Xo2.over2;
                Xo_3=Xo3.over3;
                
                
                for i=1:size(Xo_1,2)
                    for j=1:size(Xo_2,2)
                        if (max(Xo_1{i}(1,:),Xo_2{j}(1,:))==X_Over{k}(1,1:end))&(min(Xo_1{i}(2,:),Xo_2{j}(2,:))==X_Over{k}(2,1:end))
                            g9=Xo_1{i};
                            g10=Xo_2{j};
                            [dp1,dp2]=OCN_ActivationFun(Ah,g9,g10,X_Over{k}(1,:),X_Over{k}(2,:),gamma,n_fea);
                            dp1_c1=dp1;
                            dp2_c2=dp2;
                        end
                    end
                end
                
                for i=1:size(Xo_1,2)
                    for j=1:size(Xo_3,2)
                        if (max(Xo_1{i}(1,:),Xo_3{j}(1,:))==X_Over{k}(1,1:end))&(min(Xo_1{i}(2,:),Xo_3{j}(2,:))==X_Over{k}(2,1:end))
                            g9=Xo_1{i};
                            g10=Xo_3{j};
                            [dp1,dp2]=OCN_ActivationFun(Ah,g9,g10,X_Over{k}(1,:),X_Over{k}(2,:),gamma,n_fea);
                            dp1_c1=dp1;
                            dp3_c3=dp2;
                        end
                        
                    end
                end
                
                for i=1:size(Xo_2,2)
                    for j=1:size(Xo_3,2)
                        if (max(Xo_2{i}(1,:),Xo_3{j}(1,:))==X_Over{k}(1,1:end))&(min(Xo_2{i}(2,:),Xo_3{j}(2,:))==X_Over{k}(2,1:end))
                            g9=Xo_2{i};
                            g10=Xo_3{j};
                            [dp1,dp2]=OCN_ActivationFun(Ah,g9,g10,X_Over{k}(1,:),X_Over{k}(2,:),gamma,n_fea);
                            dp2_c2=dp1;
                            dp3_c3=dp2;
                        end
                    end
                end
            end
        end
        disp(['dp1_c1=',num2str(dp1_c1)]);
        disp(['dp2_c2=',num2str(dp2_c2)]);
        disp(['dp3_c3=',num2str(dp3_c3)]);
        
        
        
        %% Result
        
        c1=max(bj); %eq-5 pg-45 for class-1
        c2=max(bk); %eq-5 pg-45 for class-2
        c3=max(bl); %eq-5 pg-45 for class-3
        Co1=dp1_c1;
        Co2=dp2_c2;
        Co3=dp3_c3;
        Cc1=e_ccn_c1;
        Cc2=e_ccn_c2;
        Cc3=e_ccn_c3;
        o1=min(Co1,Cc1);    %eq-6 pg-45 for class-1
        o2=min(Co2,Cc2);    %eq-6 pg-45 for class-2
        o3=min(Co3,Cc3);    %eq-6 pg-45 for class-3
        disp(['o1=',num2str(o1)]);
        disp(['o2=',num2str(o2)]);
        disp(['o3=',num2str(o3)]);
        mu1=c1+o1;  %eq-4 pg-45 for class-1
        mu2=c2+o2;  %eq-4 pg-45 for class-1
        mu3=c3+o3;  %eq-4 pg-45 for class-1
        disp(['mu1=',num2str(mu1)]);
        disp(['mu2=',num2str(mu2)]);
        disp(['mu3=',num2str(mu3)]);
        if(mu1>mu2)&&(mu1>mu3)
            fprintf('class-1\n');
            count1(i1)=count1(i1)+1;
            i1=i1+1;
        elseif(mu2>mu1)&&(mu2>mu3)
            fprintf('class-2\n');
            count2(i1)=count2(i1)+1;
            i1=i1+1;
        elseif(mu3>mu1)&&(mu3>mu2)
            fprintf('class-3\n');
            count3(i1)=count3(i1)+1;
            i1=i1+1;
        else
            fprintf('NoN\n');
            count4(i1)=count4(i1)+1;
            i1=i1+1;
        end
        fprintf('\n');
    end
    
    
    %% % Error
    
    count1;
    count2;
    count3;
    count4;
    error_1=0;
    error_2=0;
    error_3=0;
    for i=1:n_test1
        if(count1(i)==0)
            error_1=error_1+1;
        end
    end
    for i=(n_test1+1):(n_test1+n_test2)
        if(count2(i)==0)
            error_2=error_2+1;
        end
    end
    for i=(n_test1+n_test2+1):n_test
        if(count3(i)==0)
            error_3=error_3+1;
        end
    end
    error=error_1+error_2+error_3;
    error_rst(th)=(error/(n_test))*100;
    
    
    %% No of Hyperbox Created
    
    if(X_Cont{1}==zeros(2,nFeature))
        N_cont_hyp=0;
    else
        N_cont_hyp=size(X_Cont,2);
    end
    
    if(X_Over{1}==zeros(2,nFeature))
        N_over_hyp=0;
    else
        N_over_hyp=size(X_Over,2);
    end
    
    
    N_hyp_c1=size(x_hyp,2);
    N_hyp_c2=size(y_hyp,2);
    N_hyp_c3=size(z_hyp,2);
    
    
    N_hyp(th)=N_hyp_c1+N_hyp_c2+N_hyp_c3+N_cont_hyp+N_over_hyp;  % total no of hyperbox created in hidden layer
    
    
end
end 
toc
%% Display Figure and Result

disp(error_rst);
plot(theta,error_rst ,':b o');  %performance plot
grid on
xlabel('Expansion Coefficient');
ylabel('% Error');
title('Performance Plot')
% if(strcmp(str,'Train'))
%     title('Training/learning plot');
% else
%     title('Performance plot');
% end

figure;
plot(theta,N_hyp); %hyperbox plot
grid on
xlabel('Expansion Coefficient');
ylabel('No 0f hyperbox created in hidden layer');
title('Hyperbox plot');

% hold on


