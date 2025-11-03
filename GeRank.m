function [fstat,gstat] = GeRank(model,genelist,met_id,max_loop)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

initCobraToolbox;
id_biomass=find(model.c);
id_oxygen=find(strcmp(model.rxns,"EX_o2_e"));
id_carbon=find(strcmp(model.rxns,"EX_glc__D_e"));
model.lb(id_biomass)=0.05;
model.lb(id_carbon)=-15;
model.lb(id_oxygen)=-20;
model.ub(id_carbon)=0;
model.ub(id_oxygen)=0;
opt=optimizeCbModel(model);
TMGR=opt.f;
numMultiStrat=10;
num_gene=numel(genelist);
num_loop=power(2,num_gene);
stat=zeros(num_loop,4);
kogene=cell(num_loop,1);
glist=zeros(1,num_gene);
for gi=1:num_gene
    glist(gi)=find(strcmp(model.genes, genelist(gi)));
end

% glc__D_e 1203 o2_e 1250
parfor k=1:num_loop
    [new_model,id_target,TMPR] = introExchange(model,id_biomass,[id_carbon,id_oxygen],met_id);
    tStart=tic;
    if TMPR>0
        [x_target,alpha,knockout] = ratioSelect(new_model,id_biomass,id_target,TMGR,max_loop,20*TMPR/max_loop,numMultiStrat,glist,k);
        stat(k,:)=[x_target,alpha,TMPR,toc(tStart)];
        kogene{k,1}=knockout;
    else
        x_target=0;
        alpha=0;
        stat(k,:)=[x_target,alpha,TMPR,toc(tStart)];
    end
end

if any(stat(:,1)>0.001)
    row=find(stat(:,1)==max(stat(:,1)));
    row=row(1);
    fstat=stat(row,:);
    gstat=kogene{row,1};
end

end

