function [signal_bipolar] = make_bipolar_signal(signal,ops)
%MAKE_BIPOLAR_SIGNAL Summary of this function goes here
%   Detailed explanation goes here

signal_bipolar= double([]);
%signal_bipolar_ca= double([]);
pbar=ProgressBar(size(ops.bip_ch_id_valid,1));
for bipolar_id=1:size(ops.bip_ch_id_valid,1)
    bipol_ch_1=ops.bip_ch_id_valid(bipolar_id,1);
    bipol_ch_2=ops.bip_ch_id_valid(bipolar_id,2);
    bip_ch_name1=ops.bip_ch_label_valid{bipolar_id,1};
    bip_ch_name2=ops.bip_ch_label_valid{bipolar_id,2};
    
    A =cell2mat(extract(bip_ch_name1,lettersPattern));
    B = cell2mat(extract(bip_ch_name2,lettersPattern));
    assert(all(A==B));
    B =str2num(cell2mat(extract(bip_ch_name2,digitsPattern)));
    A =str2num(cell2mat(extract(bip_ch_name1,digitsPattern)));
    assert(A-B==1)
    signal_bipolar(:,bipolar_id)=signal(:,bipol_ch_1)-signal(:,bipol_ch_2);
    pbar.step([],[],[]);
end


end

