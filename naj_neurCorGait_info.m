function procInfo = naj_neurCorGait_info(sii,file_out);

if isfile(file_out)
    load(file_out);
else
    %% store processing results
    procInfo = [];
    for si = sii
        procInfo(si).ID  = sprintf('sub-%03d', si);
    end
    save(file_out, 'procInfo');
end