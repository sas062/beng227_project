opts = struct('format','doc','showCode',true,'evalCode',false);

publish('main.m', opts)
publish('glucagon_pde_rhs.m', opts)
publish('glucagon_secretion.m', opts)
publish('insulin_secretion.m', opts)
publish('ren_model.m', opts)
publish('report_plots.m', opts)
publish('rhs_islet.m', opts)
publish('assign_cell_type.m', opts)
publish('circle_packing.m', opts)
publish('get_NN.m', opts)
