%% clear command window

clc;
clear all;
close all;

%% select input
%[fname, fpath] = uigetfile('*.mat', 'choose dataset'); %'.mat'
%if ~isequal(fname, 0) || ~isequal(fpath, 0)
%    full_path_file = fullfile(fpath, fname);
 %   load (full_path_file);
    data=load('Samples_CIELab.txt');
    samples = data;
    x = samples(:,1:end-1);          %xdata
    x_label= samples(:,end);
    %NF=max(max(xdata));             %normalization
    %x = xdata./NF;
    
    %% hyperbox Creation & Expansion
    theta = 0.01:0.02:0.2;
    l_theta=length(theta);
    gamma=4;
    
    start_time = tic;
    for rep = 1 : 10
        
        [Train, Test] = crossvalind('HoldOut', x_label, 0.1);
        TrainingSample= x(Train,:);
        TrainingLabel = x_label(Train,:);
        TestingSample = x(Test,:);
        TestingLabel  = x_label(Test,:);
        
        % ------------ modification  ------------------ 28/dec/2016
        [class_count class_label] = hist(TrainingLabel,unique(TrainingLabel));
        n_class = max(class_label);
        
        for i=1:n_class
            input1=TrainingSample((TrainingLabel== i),:);
            input(i).in = input1;
            nTrainData(i)=size(input1,1);
            
            test1=TestingSample((TestingLabel == i),:);
            test(i).in= test1;
        end
        
        nFeature=size(x,2);
        
        for th=1:l_theta
            for i = 1:n_class
                V(i).temp = zeros(nTrainData(i),nFeature);
                W(i).temp = zeros(nTrainData(i),nFeature);
                V_new(i).temp = zeros(nTrainData(i),nFeature);
                W_new(i).temp = zeros(nTrainData(i),nFeature);
            end
            
            for j=1:n_class
                V(j).temp(1,:)=input(j).in(1,:);
                W(j).temp(1,:)=input(j).in(1,:);
                HB{j}=[V(j).temp(1,:);W(j).temp(1,:)];        %initial hyperboxes one for each class
            end
            
            count=1;
            %ah=zeros(size(TrainingSample,1),nFeature);
            for  m=1:size(TrainingSample,1)
                V(m).temp(1,:)= TrainingSample(m,:);
                W(m).temp(1,:)= TrainingSample(m,:);
                HB_temp{m}=[V(m).temp(1,:); W(m).temp(1,:)];  %temp hyperbox to read input data
                %ah(m,:) = TrainingSample(m,:);
                
                for p=1:n_class
                    if(TrainingLabel(m)==class_label(p))
                        q = HB_expansion(HB{p},HB_temp{m});     %ah(m)
                        
                        if(nFeature*theta(th)>= q)&&(isolation1(HB(p),HB_temp(m),nFeature)==1)
                            %for i=1:n_class
                            %if( (i+1<= n_class) && (isolation1(HB(i),HB(i+1),nFeature)==1))
                            V_new(m).temp(1,:)= min(HB{p},[],1);
                            W_new(m).temp(1,:)= max(HB_temp{m},[],1);                %ah(m);
                            class(p).HB{count}=[V_new(m).temp(1,:); W_new(m).temp(1,:)];
                            
                        else
                            count = count+1;
                            V_new(m).temp(1,:)= V(m).temp(1,:); %;ah;
                            W_new(m).temp(1,:)= W(m).temp(1,:); %ah;
                            class(p).HB{count}=[V_new(m).temp(1,:);W_new(m).temp(1,:)];
                        end
                        HB{p}= class(p).HB{count};
                    end
                end
            end
            
            for HBC=1:n_class
                CLASS_HB =class(HBC).HB;
                emptycells = cellfun('isempty',CLASS_HB);
                CLASS_HB(emptycells)=[];
                CLASS{HBC} = CLASS_HB;                 %sum(~cellfun('isempty',class(HBC).HB));
                NHBCS(th).HB(HBC) = length(CLASS_HB);  %class(HBC).HB(~cellfun('isempty',class(HBC).HB));
            end
            for THB = 1:n_class
                Total_HBs(THB)= NHBCS(th).HB(THB);
            end
            Total_Hyperboxes(th)= sum(Total_HBs);
            
            % containment & overlap compensate neuron
            Vc = zeros(1,nFeature);
            Wc = zeros(1,nFeature);
            Vo = zeros(1,nFeature);
            Wo = zeros(1,nFeature);
            
            a=0;
            a1=0;
            
            O1=[];%zeros(2,nFeature);
            O2=[];%zeros(2,nFeature);
            
            Xc=[];
            Xo=[];
            
            for pp=1:n_class
                if pp<n_class
                    count1= NHBCS(th).HB(pp);
                    count2= NHBCS(th).HB(pp+1);
                    for q1=1:count1
                        for q2=1:count2
                            if(isolation_test(CLASS{pp}{q1},CLASS{pp+1}{q2})==1)
                                %                         Xc=[];
                                %                         Xo=[];
                                fprintf('No Overlap Or Containment\n');
                            elseif(Cont_Comp_Neuron(CLASS{pp}{q1},CLASS{pp+1}{q2})== 1)
                                a=a+1;
                                C1.HB{a}= CLASS{pp}{q1};
                                C2.HB{a}= CLASS{pp+1}{q2};
                                contain_class(a).HB =[pp pp+1];
                                for d=1:nFeature
                                    %                             Vc(d)= max(O1(1,d),O2(1,d));
                                    %                             Wc(d)= min(O1(2,d),O2(2,d));
                                    Vc(d)= max(C1.HB{1,a}(1,d),C2.HB{1,a}(1,d));
                                    Wc(d)= min(C1.HB{1,a}(2,d),C2.HB{1,a}(2,d));
                                end
                                Xc.HB{a}=[Vc;Wc];
                                
                                % Vc = max([max(CLASS{pp+1}{q2})],[max(CLASS{pp}{q1})]);
                                % Wc = min([min(CLASS{pp+1}{q2})],[min(CLASS{pp}{q1})]);
                                
                            else
                                a1=a1+1;
                                O1.HB{a1}=CLASS{pp}{q1};
                                O2.HB{a1}=CLASS{pp+1}{q2};
                                overlap_class(a1).HB=[pp pp+1];
                                for d=1:nFeature
                                    Vo(d)= max(O1.HB{1,a1}(1,d),O2.HB{1,a1}(1,d));
                                    Wo(d)= min(O1.HB{1,a1}(2,d),O2.HB{1,a1}(2,d));
                                end
                                
                                Xo.HB{a1}=[Vo;Wo];
                                
                                %Vo= max([max(CLASS{pp+1}{q2})],[max(CLASS{pp}{q1})])
                                %Wo= min([min(CLASS{pp+1}{q2})],[min(CLASS{pp}{q1})])
                                
                                
                            end
                        end
                    end
                end
            end
            
            nCCN=length(Xc);  %size(Xc.HB,2);
            nOCN=length(Xo);  %size(Xo.HB,2);
            
            
            memb_TempValue = zeros(size(x,1),n_class);   %TestingSample %TrainingSample %x
            for r=1:size(x,1)                            %x,1) %TestingSample,1
                Ah= x(r,:);                            %x(r,:);%TrainingSample(r,:); TestingSample(r,:);
                for i=1:n_class
                    for j=1:NHBCS(th).HB(i)
                        Memb_value(i).HB(j)= MembershipFunc(Ah,CLASS{i}{j}(1,:),CLASS{i}{j}(2,:),gamma);
                        cl(i) = max(Memb_value(i).HB);
                    end
                end
                memb_TempValue(r,:)=cl;
                
                ccn_value=zeros(1,n_class);
                ej_ccn=0;
                for k=1:nCCN
                    cou=0;
                    for d=1:nFeature
                        if(Ah(1,d)<=Xc.HB{1,k}(2,d) & Ah(1,d)>=Xc.HB{1,k}(1,d))
                            cou=cou+1;
                        end
                    end
                    if(cou==nFeature)
                        ej_ccn(k) = CCN_ActivationFun(Ah,Xc.HB{1,k}(1,:),Xc.HB{1,k}(2,:),gamma);
                        ccnclass1(k)= contain_class(k).HB(1);
                        %ccnclass2(k)= contain_class(k).HB(2);
                        ccn_value(ccnclass1)= ej_ccn;
                    end
                    
                end
                Ovm_value=zeros(1,n_class);
                dp1=0;
                dp2=0;
                for k=1:nOCN
                    cou=0;
                    for d=1:nFeature
                        if(Ah(1,d)<=Xo.HB{1,k}(2,d) | Ah(1,d)>=Xo.HB{1,k}(1,d)) % before &
                            cou=cou+1;
                        end
                    end
                    if(cou==nFeature)
                        [dp1,dp2]= OCN_ActivationFun(Ah,O1.HB{k},O2.HB{k},Xo.HB{1,k}(1,:),Xo.HB{1,k}(2,:),gamma,nFeature);
                        overlapclass1(k)=overlap_class(k).HB(1);
                        overlapclass2(k)=overlap_class(k).HB(2);
                        Ovm_value(overlapclass1) = dp1;
                        Ovm_value(overlapclass1) = dp2;
                    end
                    
                    %       Ovm_value(overlapclass1) = dp1;
                    %       Ovm_value(overlapclass1) = dp2;
                end
                memb_TempValue(r,:)= memb_TempValue(r,:)+Ovm_value+ccn_value;
            end
            
            for Rk=1:size(x,1)                                      %size(TestingSample,1) %x,TrainingSample
                [max_vlaue(Rk) ix(Rk)]=max(memb_TempValue(Rk,:));
            end
            
            Predicted_labels=ix;
            [C,order] = confusionmat(x_label',Predicted_labels);     %TestingLabel% TrainingLabel %x_label
            Rate_Recog(th) = (trace(C)/size(x_label,1))*100;         %TestingSample %x %TrainingLabel
            Error_Rate(th) = 100-Rate_Recog(th);
        end
        
        figure(rep);
        plot(theta,Error_Rate,':b *');
        grid on
        xlabel('Expansion Coefficient');
        ylabel('Error Rate');
        title('Performance Plot');
        drawnow
        
        figure(rep + 10);
        plot(theta,Total_Hyperboxes,':r o')
        grid on
        xlabel('Expansion Coefficient');
        ylabel('No 0f hyperbox created in hidden layer');
        title('Hyperbox plot For Iris Data');
        drawnow
    end
    
    % total time consumed
    time_consumed = toc(start_time);
    fprintf('Time consumed is %.2f second\n', time_consumed);
    
%end