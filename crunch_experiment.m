
experiment_name = 'MITLangloc';
subjects = {'AMC096'};

for i=1:length(subjects)
    crunch_subject_ALBANY(subjects{i},experiment_name)
end

%NOTE: AMC088 has MITLangloc but not MITNLengthSentences
%'AMC082', 'AMC086','AMC091','AMC092',,'AMC097','AMC099'