a = dir('*.mat');
SCE=cell(size(a,1),1);
for k=1:size(a,1)
    load(sprintf('%s',a(k).name));
    grpid = strrep(a(k).name,'.mat','');
    %sce.c_batch_id = string(repmat(grpid, sce.NumCells,1));
    %save(sprintf('mat_files/%s',a(k).name),'sce','-v7.3');
    SCE{k}=sce;
end
sce=sc_mergesces(SCE);


sce = sce.qcfilter(1000, 0.15, 10);
alternative = false;
if alternative
    sce = umap_sc_wrapper(sce,2,false,2.0);
    %sce = leiden_clustering_ann(sce,2.0,false,'mouse');
    sce = sce.assigncelltype('mouse',false);
else
    sce = sce.embedcells('umap2d',true,false,2);
    sce = sce.clustercells(50,'kmeans',true);
    sce = sce.assigncelltype('mouse',false);
end

save('../nr4a1_wholebody_ko_all_samples.mat','sce','-v7.3');


