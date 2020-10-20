experiment_name = 'MITNLengthSentences';
subjects = {'AMC096','AMC097','AMC099'}

for i=1:length(subjects)
    crunch_subject_ALBANY(subjects{i},experiment_name)
end
    
%'AMC082','AMC083', 'AMC086','AMC091','AMC092',