function PrepDataForR()
% This function generates three datasets for each condition. Each dataset
% contains all the subjects' data. Type of data that will be collected
% together should be set below (variable 'whichData')

% listOfSubjects = {'ag', 'az' , 'cs', 'gc', 'hr', 'ls', 'mp', 'sb', 'tc', 'vs'};
listOfSubjects = {'ag', 'az' , 'cs', 'ls', 'hr', 'mp', 'sb', 'tc','vs'};

numSubjects = length(listOfSubjects);
DataFolder = 'D:\Analysis\Behavioral-Normal-Subject\Final Clean Data';
SaveFolder = 'D:\Analysis\Behavioral-Normal-Subject\For R';
whichData = 'SVC_prctl';

currentFolder = pwd;
cd(DataFolder);

for subcount = 1:numSubjects
    load([whichData,'_',listOfSubjects{subcount},'.mat']);
    eval(['DATA_2(subcount,:,:) = ',whichData,'_2;']);
    eval(['DATA_6(subcount,:,:) = ',whichData,'_6;']);
    eval(['DATA_20(subcount,:,:) = ',whichData,'_20;']);
end

cd(currentFolder);
save([SaveFolder,'\DATA_2.mat'],'DATA_2');
save([SaveFolder,'\DATA_6.mat'],'DATA_6');
save([SaveFolder,'\DATA_20.mat'],'DATA_20');
end