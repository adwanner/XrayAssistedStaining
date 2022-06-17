function my_save(filename,varargin)
    NInputs=numel(varargin)/2;
    varlist=cell(1,NInputs);
    for iinput=1:NInputs
        str=[varargin{(iinput-1)*2+1} ' = varargin{' num2str(iinput*2) '};' ];
        eval(str);
        varlist{iinput} = varargin{(iinput-1)*2+1};
    end
%     save(filename,varlist{:},'-v7.3');
    save(filename,varlist{:});
end