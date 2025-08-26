D=readtable('../2023_nr4a1_colon_meta_info.xlsx');
id = D.sampleID;
grpid = D.grpid;

for k = 1:length(id)
    if exist(string(grpid{k})+".mat","file"), continue; end
    
    f=fullfile('../h5_files_pools_1_2',num2str(id(k)),'sample_filtered_feature_bc_matrix.h5');
    %exist(f,"file")
    [X, g, b] = sc_read10xh5file(f);
    fprintf('Raw sample %s has ncells %d and ngenes %d \n', grpid{k}, size(X,2), size(X,1));
    sce=SingleCellExperiment(X,g);
    sce.c_cell_id = b;
    sce.c_batch_id = string(repmat(grpid{k},sce.NumCells,1));
    % Light QC
    %sce = onestepraw2anno(sce,'mouse');
    % Light QC
    sce = sce.qcfilter(500, 0.20, 15);
    % Stringent QC
    %sce = sce.qcfilter(1000, 0.15, 10);
    % Optional for individual visualization
    %sce = sce.embedcells('umap2d',true);
    %sce = sce.clustercells(30, 'kmeans', true);
    %sce = sce.assigncelltype('mouse', false);
    fprintf('QC sample %s has ncells %d and ngenes %d \n', grpid{k}, size(sce.X,2), size(sce.X,1));

    save(grpid{k},'sce','-v7.3');
end


