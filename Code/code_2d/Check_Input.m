function  path = Check_Input(InP, total_freqs)
% The function is to check the validity of inputs of Inp. It is necessary
% to avoid the wrong typewriting.
%=========================check input========================
check_value=0;

if InP.new_res>=100 || InP.new_res<0.5
    check_value=1;
end

% % if InP.SNR>1000
% %     check_value=1;
% % end

if InP.linear_method>10 || InP.linear_method<=0
    check_value=1;
end

if InP.model_phantom<1
    check_value=1;
end
% output
if check_value
    error('Please check the value of input!')
end
%===================check if self-defined frequencies exist in the original frequencies
load(['..\data\model' num2str(InP.model_phantom) '\all_material.mat']);
temp = cell2mat(total_freqs);
loc = find(ismembertol(temp,frequency)==0);
if ~isempty(loc)
    error(['The input of ' num2str(temp(loc)/1e9) ' GHz cannot be found in original frequencies!']);
end
% =================check and make the new folder for new test===========
path=['..\result\' InP.test_name '\'];
if ~isdir(path)
    mkdir(path);
else
    disp('the folder is existed!');
    return;
end

end
