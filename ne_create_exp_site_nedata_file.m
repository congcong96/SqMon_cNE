function exp_site_nedata = ne_create_exp_site_nedata_file(nefile)

% updated 7/6/17 by JS to include stimlength and probetype from filename

    clear('exp_site_nedata','nedata');
    load(nefile);
    
    if exist('exp_site_nedata','var') && sum(isfield(exp_site_nedata,{...
            'nedata', 'exp','site', 'stim','df','fs','probetype','depth'}))==8 ...
            && length(fieldnames(exp_site_nedata)) == 8
 
        fprintf('\nexp_site_nedata for %s already exists!\n',nefile);
     
    else
                fprintf('\nProcessing exp_site_nedata for %s...\n',nefile);
                exp_site_nedata = [];
                
                if exist('nedata','var')
                    exp_site_nedata.nedata = nedata;
                end

                [~,nefile,~] = fileparts(nefile);
                exp_site_nedata.exp = regexp(nefile,'^\d{6}_\d{6}','match','once');
                exp_site_nedata.site = str2double(regexp(nefile,'(?<=(-site))\d{1,2}(?=(-))','match','once'));
                exp_site_nedata.stim = regexp(nefile,'\w+(?=(-\d+min))','match','once');
                exp_site_nedata.df = str2double(regexp(nefile,'(?<=(ne-))\d+(?=(dft))','match','once'));
                exp_site_nedata.fs = str2double(regexp(nefile,'(?<=(-fs))\d+(?=(-))','match','once'));
                exp_site_nedata.probetype = regexp(nefile, '(?<=(min-))\w+(?=(-fs))', 'match', 'once');
                exp_site_nedata.stimlength = str2double(regexp(nefile, '(?<=(-))\d{1,3}(?=(min-))' ,'match','once'));
                exp_site_nedata.depth = str2double(regexp(nefile,'(?<=(-))\d+(?=(um-))','match','once'));
                
                % for old file names without time and probe
                if isempty(exp_site_nedata.stim)
                    exp_site_nedata.stim = regexp(nefile,'(?<=(db-))\w+(?=(-fs))','match','once');
                end
                
                if isempty(exp_site_nedata.probetype)
                    exp_site_nedata.probetype = 'a1x32poly3';
                end
                
                if isnan(exp_site_nedata.stimlength)
                    exp_site_nedata.stimlength = 10;
                end
               
    end
end


