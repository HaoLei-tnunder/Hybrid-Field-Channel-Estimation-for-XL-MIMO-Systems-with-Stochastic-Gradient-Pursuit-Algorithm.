
N = 100:100:300; 
N2 = 2 ;
M = 2;
len = length(N);
L = 16 ;

fc = 30e9; % carrier frequency
c = 3e8;
lambda_c = c/fc; % wavelength
d = lambda_c / 2; % antenna space

HF_OMP = zeros(len,1);
HF_SGP = zeros(len,1);
HF_OMP_1 = zeros(len,1);
HF_SGP_1 = zeros(len,1);


for i = 1 : len

    % the near-field polar-domain transform matrix [5]
    Rmin=10;
    Rmax=80;
    eta = 2.5; 
    [Polarcodebook, ~,~, ~] = QuaCode(N(i), d, lambda_c, eta, Rmin, Rmax);
    % [Un, label, dict_cell, label_cell] = QuaCode(N, d, lambda_c, eta, Rmin, Rmax);
    S = size(Polarcodebook, 2);


    %%  HF-OMP

    HF_OMP(i) = L * N(i)*N2 * (  N(i)*N2 + S *M   )  + N(i)*N2* ( N(i) +N2   )  + N(i) * M * ( N2 + S)  ;

    for l = 1 : L/2
        HF_OMP(i) = HF_OMP(i) +  2*  (  1/2* l^3 + 3/2*l^2 + 2*N(i)*N2*l^2 + N(i)*N2*l  ) ;
    end


    %%  HF-SGP

    HF_SGP(i) =  L * N(i)*N2 * (N(i)*N2 + S *M)  + N(i)*N2 *( N(i) +N2   )  + N(i) * M * ( N2 + S)  ;

    for l = 1 : L/2
        HF_SGP(i) = HF_SGP(i) +  2*  ( N(i)*N2*(l*2+1)  ) ;
    end


    %%  HF-OMP-1

    for l = 1 : L
        HF_OMP_1(i) = HF_OMP_1(i) +  1/2*l^3 + 3/2*l^2 + 2*N(i)*N2*l^2 + N(i)*N2*l ;
    end

    for j = 0 : (L-1)
        for l = 1 : (L-j)
            HF_OMP_1(i) = HF_OMP_1(i) +  1/2*l^3 + 3/2*l^2 + 2*N(i)*N2*l^2 + N(i)*N2*l  ;
        end
    end

    HF_OMP_1(i) = HF_OMP_1(i) +  L* N(i)^2 * N2^2 + 1/2*L*(1+L)*N(i)*N2*S*M +  L* N(i)^2 * N2^2 +  1/2*L*(1+L)*(   N(i)*N2*( N(i)*N2+  S*M     )        )      + N(i)*N2 *( N(i) +N2   )  + N(i) * M * ( N2 + S) ; 


    %%  HF-SGP-1    

    for l = 1 : L
        HF_SGP_1(i) = HF_SGP_1(i) +  ( N(i)*N2*(l*2+1)  )  ;
    end

    for j = 0 : (L-1)
        for l = 1 : (L-j)
            HF_SGP_1(i) = HF_SGP_1(i) +  ( N(i)*N2*(l*2+1)  )  ;
        end
    end

    HF_SGP_1(i) = HF_SGP_1(i) +  L* N(i)^2 * N2^2 + 1/2*L*(1+L)*N(i)*N2*S*M +  L* N(i)^2 * N2^2 +  1/2*L*(1+L)*(   N(i)*N2*( N(i)*N2+  S*M     )        )      + N(i)*N2 *( N(i) +N2   )  + N(i) * M * ( N2 + S) ; 


%     %%  HF-OMP
% 
%     HF_OMP(i) = (L+1)*(N(i)+S)*N(i) ;
%     for l = 1 : L/2
%         HF_OMP(i) = HF_OMP(i) + l^3 + 3*l^2 + 4*N(i)*l^2 + 2*N(i)*l ;
%     end
% 
%     %%  HF-SGP
% 
%     HF_SGP(i) = (L+1)*N(i)^2+(L+1)*S*N(i) + (2+L/2)*L*N(i);
% 
%     %%  HF-OMP-1
% 
%     for l = 1 : L
%         HF_OMP_1(i) = HF_OMP_1(i) + 2*N(i)^2 + 1/2*l^3 + 3/2*l^2 + 2*N(i)*l^2 + N(i)*l ;
%     end
% 
%     for j = 0 : (L-1)
%         for l = 1 : (L-j)
%             HF_OMP_1(i) = HF_OMP_1(i) + 2*S*N(i) + 1/2*l^3 + 3/2*l^2 + 2*N(i)*l^2 + N(i)*l + N(i)*N(i) ;
%         end
%     end
% 
%     HF_OMP_1(i) = HF_OMP_1(i) +   N(i)*N(i) + S*N(i) ; 
% 
% 
%     %%  HF-SGP-1    
% 
%     for l = 1 : L
%         HF_SGP_1(i) = HF_SGP_1(i) + 2*N(i)^2 + N(i)*(2*l+1) ;
%     end
% 
%     for j = 0 : (L-1)
%         for l = 1 : (L-j)
%             HF_SGP_1(i) = HF_SGP_1(i) + 2*S*N(i) + N(i)*(2*l+1) + N(i)*N(i) ;
%         end
%     end
% 
%     HF_SGP_1(i) = HF_SGP_1(i) +   N(i)*N(i) + S*N(i) ;     



end

% figure('color',[1,1,1]); hold on; box on; grid on;
% plot(  N,   HF_OMP               ,      '     k   -   ^ '    , 'linewidth',  1);
% plot(  N,   HF_OMP_1               ,      '     b   -  o  '    , 'linewidth',  1);
% plot(  N,   HF_SGP   ,      '     r  --    p'    ,  'linewidth' , 1);
% plot(  N,   HF_SGP_1   ,      '     r  --   d '    ,  'linewidth' , 1);
% % xlim([128 512])
% xlabel('Number of antennas at the BS','Interpreter','Latex');
% ylabel('Number of complex multiplications','Interpreter','Latex');
% legend(  'Hybrid-field OMP (with $\gamma$)',  'Hybrid-field OMP (without $\gamma$)','Hybrid-field SGP (with $\gamma$)',  'Hybrid-field SGP (without $\gamma$)', 'Interpreter','Latex');
% 





%%

figure('color',[1,1,1]); box on; grid on;

Y=[     HF_SGP    HF_OMP    HF_SGP_1   HF_OMP_1    ];

% %准备数据
% Y=[70,75,80,85;
%      80,85,90,95;
%      90,95,100,105;
%     100,105,110,115];
X=1:3;
 %画出4组柱状图，宽度1
h=bar(X,Y,1);      
 %修改横坐标名称、字体
set(gca,'XTickLabel',{'100','200','300'},'FontSize',10,'FontName','Times New Roman');
% 设置柱子颜色,颜色为RGB三原色，每个值在0~1之间即可
set(h(1),'FaceColor',[30,150,252]/255)     
set(h(2),'FaceColor',[162,214,249]/255)    
set(h(3),'FaceColor',[252,243,0]/255)    
set(h(4),'FaceColor',[255,198,0]/255)   
set(h(4),'FaceColor',[255,198,100]/255) 
%修改x,y轴标签
xlabel('Number of antennas at the BS','Interpreter','Latex');
ylabel('Number of complex multiplications','Interpreter','Latex');
%修改图例
legend('Hybrid-field SGP (with $\gamma$)',    'Hybrid-field OMP (with $\gamma$) [22]',        'Hybrid-field SGP (without $\gamma$)',    'Hybrid-field OMP (without $\gamma$) [26]',   'Interpreter','Latex');


