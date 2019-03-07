Xorgi= xlsread('lll.xlsx','Sheet2');
Xorg= abs(Xorgi(:,1:6));
aa=ones(size(Xorg,1),6).*repmat(max(Xorg),size(Xorg,1),1);
Xorg = Xorg./aa;

Xorg = [Xorg Xorgi(:,14)]; % Data for four classes

Xla= Xorg(Xorg(:,7)==1,1:6); % class 1 data (la tree species)
Xar= Xorg(Xorg(:,7)==2,1:6); % class 2 data (ar tree species)
Xab= Xorg(Xorg(:,7)==3,1:6); % class 3 data (ab tree species)
Xpc= Xorg(Xorg(:,7)==4,1:6); % class 4 data (pc tree species)
% 

%op = [mean(Xla(:,:)); mean(Xar(:,:)); mean(Xab(:,:)); mean(Xpc(:,:))];

%[rho, pval] = corr(op, 'rows','pairwise');

x = [Xla;Xar;Xab;Xpc]; x = x(:); 
g1 = [repmat({'la'},size(Xla)); repmat({'ar'},size(Xar)); repmat({'ab'},size(Xab)); repmat({'pc'},size(Xpc))];  g1 = g1(:); 
g2 = [repmat({'EGF1'},80,1); repmat({'EGF2'},80,1); repmat({'EGF3'},80,1); repmat({'EGF4'},80,1) ; repmat({'EGF5'},80,1) ; repmat({'EGF6'},80,1)];  
g2 = g2(:); 

h=boxplot(x, {g2,g1}, 'colorgroup',g1, 'factorgap',5, 'factorseparator',1);
set(gca,'xtickmode','auto','xticklabelmode','auto')
xlabel('Internal Geometric Features', 'FontSize', 16);
ylabel('Feature value', 'FontSize', 18);

%set(findobj(get(h(1), 'parent'), 'type', 'text'), 'fontsize', 14,'Rotation', 0);

%xtix = {'la','ar','ab','pc','','la','ar','ab','pc','','la','ar','ab','pc','','la','ar','ab','pc','','la','ar','ab','pc','','la','ar','ab','pc',''};   % Your labels
%xtixloc = [1:0.9:30];      % Your label locations

%set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);

% text_h = findobj(gca, 'Type', 'text');
% tickLabelStr = {'','','IGF1','','','','IGF2','','','','IGF3','','','','IGF4','','','','IGF5','','','','IGF6',''} ; 
% for cnt = 1:length(text_h)
%     set(text_h(cnt),    'FontSize', 12,...
%                         'Rotation', 0, ...
%                         'String', tickLabelStr{length(tickLabelStr)-cnt+1}, ...
%                         'HorizontalAlignment', 'right')
% end
%     
% x = [Xla;Xar;Xab;Xpc]; x = x(:); 
% %g1 = [ones(size(Xla)); 2*ones(size(Xar)); 3*ones(size(Xab)); 4*ones(size(Xpc))];  g1 = g1(:); 
% g1 = [repmat({'la'},size(Xla)); repmat({'ar'},size(Xar)); repmat({'ab'},size(Xab)); repmat({'pc'},size(Xpc))];  g1 = g1(:); 
% g2 = [repmat([repmat({'F1'},1,1); repmat({'F2'},1,1); repmat({'F3'},1,1); repmat({'F4'},1,1) ; repmat({'F5'},1,1) ; repmat({'F6'},1,1) ; repmat({'F7'},1,1) ; repmat({'F8'},1,1) ; repmat({'F9'},1,1) ; repmat({'F10'},1,1) ; repmat({'F11'},1,1) ; repmat({'F12'},1,1)],54,1)];  g2 = g2(:); 
% %g2 = repmat(1:12,54,1); g2 = g2(:);
% 
% boxplot(x, {g2,g1}, 'colorgroup',g1, 'factorgap',5, 'factorseparator',1) 
% 


%group = [repmat({'First'}, size(Xla,1), 1); repmat({'Second'}, size(Xar,1),1); repmat({'Third'}, size(Xab,1),1); repmat({'Fourth'},  size(Xpc,1),1)];
%boxplot([Xla(:,1:2);Xar(:,1);Xab(:,1);Xpc(:,1)], group);

%  op =[];
%  for i=1:1:1%size(Xla,2)
%      op = [op Xla(:,i) Xar(:,i) Xab(:,i) Xpc(:,i)];
%  end
% 
% 
% %subplot(2,2,1);
% boxplot(Xla,'factorgap',[5 4])

% subplot(2,2,2);
% boxplot(Xar)
% 
% subplot(2,2,3);
% boxplot(Xab)
% 
% subplot(2,2,4);
% boxplot(Xpc)

% Xorg= xlsread('Feature_LiDAR.xlsx','Sheet5');
% Xorg= Xorg(:,1:12);
% Xorg= abs(Xorg);
% aa=ones(size(Xorg,1),12).*repmat(max(Xorg),size(Xorg,1),1);
% Xorg = Xorg./aa;
% subplot(2,1,2);
% boxplot(Xorg)
